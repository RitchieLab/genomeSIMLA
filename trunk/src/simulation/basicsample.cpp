//
// C++ Implementation: population
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "utility/types.h"
#include "basicsample.h"
#include "poolmanager.h"
#include "pedigreesample.h"
#include "meantable.h"
#include <iomanip>

namespace Simulation {

using namespace std;
using namespace Utility;


void Sample::ApplyMissingData(PoolManager *mgr) {
	uint sampleSize = people.size();

	if (sampleSize < 0)
		return;
	PoolManager::Iterator itr = mgr->GetIterator();

	ChromPool *ch = itr.GetNext();
	//Iterate over each of the chromosomes
	while (ch) {
		//uint chrCount = ch->GetID();

		uint locusSize = ch->GetLociCount();

		//Iterate through each locus and modify genoError% of individuals
		for (uint locus=0; locus<locusSize; locus++) { 
			uint errorCount = (uint)((float)sampleSize * missingData);

			Utility::BitSetType alreadyUsed(sampleSize, false);

			while (errorCount > 0) {
				int idx = Utility::Random::globalGenerator((int)sampleSize);
				if (!alreadyUsed[idx]) {
					people[idx]->ClearLocus(ch->GetID(), locus);
					alreadyUsed[idx]=true;
					errorCount--;
				}
			}
		}
		ch = itr.GetNext();
	}
}

	

void Sample::ApplyErrors(PoolManager *mgr) {
	ApplyGenocopyError( mgr );
	//ApplyPhenocopyError( mgr );
	ApplyMissingData( mgr );
}

/**
 * isFemale(0) is male
 */
Individual *Sample::DrawIndividual(int id, int pedID, PoolManager *pools, PenetranceModel* model, int status, bool isFemale) {
	Individual* ind = NULL;

	while (ind == NULL) {
		Individual *individual = new Individual(id, pedID, pools->GetPoolCount());
		pools->DrawIndividual(*individual, isFemale);
		if (model){
			individual->ApplyStatus( model );
			if (status == -1 || individual->GetStatus() == (uint)status)
				ind = individual;
			else
				delete individual;
		}
	}
	return ind;
}

void BasicSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model) {
	uint sampleSize = people.size();
	uint errorCount = (uint)((float)affCount * phenoError);
	uint topID = people[sampleSize - 1]->GetID();
	
	vector<int> indices;
	for (uint i=0; i<sampleSize; i++) 
		indices.push_back(i);
	random_shuffle(indices.begin(), indices.end(), Utility::Random::globalGenerator);

	BitSetType tested(sampleSize, false);

	Individual *individual=NULL;
	int n=0;

	for (uint i=0; i<errorCount; i++) {
		bool isAffected = false;
		int idx;
		while (!isAffected) {
			idx = indices[n++];
			if (!tested[idx]) {
				individual = people[idx];
				isAffected = individual->IsAffected();
				tested[idx]=true;
			}
		}

		
		assert(individual->IsAffected());
		delete individual;
		bool isFemale = Random::globalGenerator.drand() < Individual::FtoM_BirthRatio;
		individual = DrawIndividual(++topID, topID, pools, model, 0, isFemale);
assert(!individual->IsAffected());
		individual->SetStatus( true );
		people[idx] = individual;
	}

}

bool Sample::ResolveGenotypes(ChromPool *pool, size_t chrID, bool justLoadedChrom) {
	vector<Individual *>::iterator itr = people.begin();
	vector<Individual *>::iterator end = people.end();
	bool success = true;

	vector<uint> fakeLoci;

	while (itr != end) {
		if (justLoadedChrom)
			pool->ResolveGenotypes( *itr );
		else
			(*itr)->ResolveGenotypes(chrID, fakeLoci);
		itr++;
	}
	return success;
}

int Sample::WriteDataset(const char *filename, uint *gtCount) {
	ofstream file(filename);
	uint observations = WriteDataset(file, gtCount);
	file.close();	
	return observations;
}
/**
 * Apply the error according to the locus rate
 */
void Sample::ApplyGenocopyError(PoolManager *mgr) {
	uint sampleSize = people.size();

	if (sampleSize < 0)
		return;
	PoolManager::Iterator itr = mgr->GetIterator();

	ChromPool *ch = itr.GetNext();
	//Iterate over each of the chromosomes
	while (ch) {
//		uint chrCount = ch->GetID();

//		for (uint chrID=0; chrID<chrCount; chrID++) {
			uint locusSize = ch->GetLociCount();

			//Iterate through each locus and modify genoError% of individuals
			for (uint locus=0; locus<locusSize; locus++) { 
				//float genoError = ch->GetErrorRate(locus);
				uint errorCount = (uint)((float)sampleSize * genoError);
				Utility::BitSetType alreadyUsed(sampleSize, false);
	
				while (errorCount > 0) {
					uint idx = Utility::Random::globalGenerator((int)sampleSize);
					if (!alreadyUsed[idx]) {
						people[idx]->ChangeGenotype(ch->GetID(), locus, 0);
						alreadyUsed[idx]=true;
						errorCount--;
					}
				}
			}
//		}
		ch = itr.GetNext();
	}
}

void Sample::Purge() {
	uint count = people.size();
	for (uint i=0; i<count; i++)
		delete people[i];
	people.clear();
}



void BasicSample::BuildSample(PoolManager &pools, PenetranceModel *model) {
	uint affecteds 		= affCount;			//Affecteds/lower tail
	uint unaffecteds 	= unaffCount;		//unaffecteds/upper tail
	uint midrange  		= midCount;			//Mid section (only on continuous)
	affCount = 0;

	//Let's inspect the first 10 chromosomes in the pool. Something is killing the
	//Windows version...gotta figure out where
	//ChromPool *pool = pools.GetIterator().GetNext();
/*	pool->SamplePool(20, 20);
	cerr<<"--SamplePool(10, 20)\n";
*/	
	uint id=0;					///<We have to keep up with individuals separately

	bool isViable = false;
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models
	while (affecteds > 0 || unaffecteds > 0 || midrange > 0) {
		bool validIndividual=false;
		isViable=false;
		while (!validIndividual) {
			Individual *newPerson = new Individual(id, id, pools.GetPoolCount());
			pools.DrawIndividual(*newPerson);
/*
			PoolManager::Iterator itr = pools.GetIterator();
			ChromPool *pool = itr.GetNext();
#ifdef USE_XY
			pool
			AlleleSource<LocusXY> *x = pools.DrawX();

			if (Utility::Random::globalGenerator.drand() < 0.5) {
				AlleleSource<LocusXY> *x2 = pools.DrawX();
				newPerson->SetXX(x, x2);
			}
			else	{
				AlleleSource<LocusXY> *y = pools.DrawX();
				newPerson->SetXY(x, y);
			}
#endif 
			while (pool) {
				pool->ApplyChromosomeIDs(newPerson);
				//At this point, the pools may or may not be in memory
				pool->ApplyPresentGenotypes(newPerson, true);
				pool=itr.GetNext();
			}
*/
			//At this point, we have all the data associated with our person
			//Time to assign status 
			int status = newPerson->ApplyStatus(model);
			PoolManager::Iterator itr = pools.GetIterator();
			ChromPool *pool = itr.GetNext();
			if (status == 1) {
				if (affecteds > 0) {
					id++;
					affecteds--;
					people.push_back(newPerson);
					validIndividual = true;
					affCount++;
					pool = itr.GetNext();
					while (pool) {
						pool->ReserveIndividual( newPerson);
						pool=itr.GetNext();
					}
					isViable =true;
				}
			}
			else if (status == 0) {
				if (unaffecteds > 0) {
					id++;
					unaffecteds--;
					people.push_back(newPerson);
					validIndividual = true;

					pool = itr.GetNext();
					while (pool) {
						pool->ReserveIndividual( newPerson);
						pool=itr.GetNext();
					}

					isViable=true;
				}
			}
			//Midrange for continuous
			else if (status == 2) {
				if (midrange > 0) {
					id++;
					midrange--;
					people.push_back(newPerson);
					validIndividual = true;

					pool = itr.GetNext();
					while (pool) {
						pool->ReserveIndividual( newPerson);
						pool=itr.GetNext();
					}

					isViable=true;
				}

			}
			if (!isViable) {
				itr=pools.GetIterator();
				pool=itr.GetNext();
				delete newPerson;
			}
		}
	}
/*
	PoolManager::Iterator itr = pools.GetIterator();
	pool = itr.GetNext();
	while (pool) {	
		//This will call populate each individual with their genotypes
		//If the pool is asleep, it will be awakened first, then put back to sleep
		pool->ResolveGenotypes(people);
		pool = pool->GetNext();
		
	}
*/
}

/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 * @note if percentAffected is 0.0, we will not be writing the status
 */
int ContinuousSample::WriteDataset(ostream &os, uint *gtCounts) {
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WriteContinuous(os, gtCounts, affCount > 0);

	return count;
}
int ContinuousSample::WriteBinaryDataset(ostream &metadata, ostream &genotypes) {
	WriteBinaryHeader(genotypes, 1, "CC  ");
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WriteContinuous(metadata, genotypes, missingData > 0.0, affCount > 0);
	return count;
}

void ContinuousSample::AppendGenotypeCountsToReport(StatusModel::ModelLociArray loci) {
	uint indCount = people.size();
	uint modelSize = loci.size();
	uint genotypeCount = (uint)pow(3.0, (int)modelSize);

	
	if (overallAffectedCounts == NULL) {
		overallAffectedCounts = new uint[genotypeCount];
		overallUnaffectedCounts = new uint[genotypeCount];
		overallOtherCounts		= new uint[genotypeCount];
		for (uint i=0; i<genotypeCount; i++) {
			overallAffectedCounts[i]=0;
			overallUnaffectedCounts[i]=0;
			overallOtherCounts[i]=0;
		}
	}
	vector<uint> genotypes;
	vector<StatusModel::DiseaseLocus>::iterator itr;
	vector<StatusModel::DiseaseLocus>::iterator end = loci.end();

	for (uint ind = 0; ind<indCount; ind++) {
		Individual &individual = *people[ind];
		itr = loci.begin();

		genotypes.clear();

		while (itr != end)  {
			genotypes.push_back(individual.GetGenotype(itr->chromosome, itr->locusIdx));
			itr++;
		}

#ifdef DEBUG_PRODUCTION
		cout<<"AppendGenotypeCountsToReport ("<<ind<<")\t";
#endif

		uint idx = BuildGenotypeIndex(genotypes);
		if (individual.GetStatus() == 0)
			overallAffectedCounts[idx]++;
		else if (individual.GetStatus() == 1)
			overallUnaffectedCounts[idx]++;
		else 
			overallOtherCounts[idx]++;
	}	
	
	
}

void ContinuousSample::ReportGenotypeCounts(vector<StatusModel::DiseaseLocus> loci, ostream& output) {
	uint modelSize = loci.size();
	
	vector<string> labels;
	BuildGenotypeLabels(labels, modelSize);
	uint colWidth = 8;
	if (labels[0].length() > colWidth)
		colWidth = labels[0].length() + 2;

	output<<setw(12)<<"Status "<<"|";
	for (uint i=0; i<labels.size(); i++) 
		output<<setw(colWidth)<<labels[i];
	output<<"\n";
	
	output<<setw(12)<<"Left Tail "<<"|";
	for (uint i=0; i<labels.size(); i++) 
		output<<setw(colWidth)<<overallAffectedCounts[i];
	output<<"\n";
	output<<setw(12)<<"Right Tail "<<"|";
	for (uint i=0; i<labels.size(); i++) 
		output<<setw(colWidth)<<overallUnaffectedCounts[i];
	output<<"\n";
	output<<setw(12)<<"Mid Section "<<"|";
	for (uint i=0; i<labels.size(); i++) 
		output<<setw(colWidth)<<overallOtherCounts[i];
	output<<"\n";

	if (overallAffectedCounts)
		delete[] overallAffectedCounts;
	if (overallUnaffectedCounts) 
		delete[] overallUnaffectedCounts;
	if (overallOtherCounts)
		delete[] overallOtherCounts;

	overallAffectedCounts 	= NULL;
	overallUnaffectedCounts = NULL;
	overallOtherCounts		= NULL;


}

/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 * @note if percentAffected is 0.0, we will not be writing the status
 */
int BasicSample::WriteDataset(ostream &os, uint *gtCounts) {
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WriteMDR(os, gtCounts, affCount > 0);

	return count;
}


uint Sample::GetMultiplier(uint genotype, uint position, uint modelSize) {
	assert(genotype < 3);
	uint pwr=modelSize - position - 1;
	return (uint)((float)genotype * powf(3.0, (float)pwr));
}

uint Sample::BuildGenotypeIndex(vector<uint> & genotypes) {
    uint numGenos = genotypes.size();
    uint index=0;

#ifdef DEBUG_PRODUCTION
	cout<<"Not converting to simla\n";
#endif
    for(uint i=0; i<numGenos; i++){
		index += GetMultiplier( genotypes[i], i, numGenos);
    }

	return index;	
}

void Sample::PopulatePool(std::vector<Chromosome>& pool, uint chrID) {
	vector<Individual*>::iterator itr = people.begin();
	vector<Individual*>::iterator end = people.end();

	while (itr != end) {
		(*itr)->PopulateChromosomePool(chrID, pool);
		itr++;
	}
}


string Sample::BuildGenotypeLabel(uint genotype, uint position) {
	char A='A'+position;
	char a='a'+position;
	
	char *label = new char[3];
	if (genotype == 0)
		sprintf(label, "%c%c", A, A);
	else if (genotype == 1)
		sprintf(label, "%c%c", A, a);
	else if (genotype == 2)
		sprintf(label, "%c%c", a, a);
	string finalLabel = label;
	delete[] label;
	return finalLabel;
}

string Sample::BuildGenotypeLabel(uint *genotypes, uint modelSize) {
	stringstream ss;
	for (uint i=0; i<modelSize; i++) 
		ss<<BuildGenotypeLabel(genotypes[i], i);
	return ss.str();
}

void Sample::BuildGenotypeLabels(vector<string>& labels, uint modelSize) {
	uint *genotypes = new uint[modelSize+1];
	uint position = 0;

//	for (uint i=0; i<modelSize; i++) 
//		genotypes[i]=0;
	memset((void*)genotypes, 0, (modelSize + 1 )*sizeof(uint));

	labels.clear();	
	position = modelSize - 1;

	while (genotypes[0]<3) {	
		string label = BuildGenotypeLabel(genotypes, modelSize);
		
		labels.push_back(label);

		if (++genotypes[position]>2 && position > 0) {
			//Find the highest position of rollover
			while (position-- > 0 && ++genotypes[position] > 2) {}

			while (position < modelSize - 1) 
				genotypes[++position] = 0;
		}
	}	
	delete[] genotypes;
}
	
void Sample::AppendGenotypeCountsToReport(StatusModel::ModelLociArray loci) {
	uint indCount = people.size();
	uint modelSize = loci.size();
	uint genotypeCount = (uint)pow(3.0, (int)modelSize);

	if (overallAffectedCounts == NULL) {
		overallAffectedCounts = new uint[genotypeCount];
		overallUnaffectedCounts = new uint[genotypeCount];
		overallOtherCounts = new uint[genotypeCount];
		
		for (uint i=0; i<genotypeCount; i++) {
			overallAffectedCounts[i]=0;
			overallUnaffectedCounts[i]=0;
			overallOtherCounts[i] =0;
		}
				
	}
	vector<uint> genotypes;
	vector<StatusModel::DiseaseLocus>::iterator itr;
	vector<StatusModel::DiseaseLocus>::iterator end = loci.end();

	for (uint ind = 0; ind<indCount; ind++) {
		Individual &individual = *people[ind];
		itr = loci.begin();

		genotypes.clear();
		while (itr != end)  {
			uint chr = itr->chromosome;
			uint loc = itr->locusIdx;
			int genotype = individual.GetGenotype(chr, loc);
			genotypes.push_back(genotype);
			itr++;
		}

#ifdef DEBUG_PRODUCTION
		cout<<"AppendGenotypeCountsToReport ("<<ind<<")\t";
#endif

		uint idx = BuildGenotypeIndex(genotypes);
		if (individual.IsAffected())
			overallAffectedCounts[idx]++;
		else
			overallUnaffectedCounts[idx]++;
	}	
	
	
}

void Sample::ReportGenotypeCounts(vector<StatusModel::DiseaseLocus> loci, ostream& output) {
	uint modelSize = loci.size();
	
	vector<string> labels;
	BuildGenotypeLabels(labels, modelSize);
	uint colWidth = 8;
	if (labels[0].length() > colWidth)
		colWidth = labels[0].length();

	output<<setw(colWidth)<<"Status "<<" | "<<setw(colWidth)<<"Affected "<<" | "<<setw(colWidth)<<"Unaffected "<<" | "<<setw(colWidth)<<"Aff-Unaff"<<"\n";
	for (uint i=0; i<labels.size(); i++) {
		output<<setw(colWidth)<<labels[i]<<" | "<<setw(colWidth)<<overallAffectedCounts[i]<<" | "<<setw(colWidth)<<overallUnaffectedCounts[i]<<" | "<<setw(colWidth)<<(int)(overallAffectedCounts[i] - overallUnaffectedCounts[i])<<"\n";
	}
	output<<"\n";
	

	if (overallAffectedCounts)
		delete[] overallAffectedCounts;
	if (overallUnaffectedCounts) 
		delete[] overallUnaffectedCounts;
	if (overallOtherCounts)
		delete[] overallOtherCounts;

	overallAffectedCounts 	= NULL;
	overallUnaffectedCounts = NULL;
	overallOtherCounts		= NULL;

}

void Sample::WriteConfiguration(ostream& file) {
	file<<GetDetails()<<"\n";
}

string BasicSample::GetDetails() {
	stringstream ss;
	ss<<"DATASET CC "<<description<<" "<<affCount<<" "<<unaffCount<<" "<<genoError<<" "<<phenoError<<" "<<missingData;
	return ss.str();
}

Sample *Sample::ParseBinaryFile(const char *gtFilename, const char *metaFilename, uint fileVersion) {
	uint version, bpb, flags;
	uint peopleCount, chromCount;

	char *buffer = new char[5];

	ifstream *genotypes, *metadata;

	genotypes = new ifstream(gtFilename);
	if (!genotypes->is_open() ) {	
		cout<<"Unable to open file: "<<gtFilename<<". Unable to continue\n";
		abort();
	}
	if (metaFilename) 
		metadata = new ifstream(metaFilename);
	else 
		metadata = NULL;


	genotypes->read(buffer, 4);
	string fileType = buffer;
	
	genotypes->read(buffer, 2);
	string endian = buffer;
	assert(strncmp(endian.c_str(), "LE", 2)==0);

	genotypes->read((char *)&version, 4);
	assert(version == fileVersion);
	genotypes->read((char *)&bpb, 4);
	assert(bpb == boost::dynamic_bitset<>::bits_per_block);
	genotypes->read((char *)&flags, 4);
	genotypes->read((char *)&peopleCount, 4);
	genotypes->read((char *)&chromCount, 4);
	
	vector<uint> chromSizes;	
	uint count;
	for (uint i=0; i<chromCount; i++) {
		genotypes->read((char *)&count, 4);
		chromSizes.push_back(count);
	}

	float missingData = 0.0;

	if (flags || 1) {
		cout<<"Missing data flag has been set....just priming to 0.1\n";
		missingData = 0.1;
	}

	Sample *theSample = NULL; 
	if (strcmp(fileType.c_str(), "PED ") == 0) {
		cout<<"Loading Binary Pedigree Data\n";
		theSample = new PedigreeSample(0.0, 0.0, missingData, "frombin");
		assert(metadata != NULL);
	} else if (strcmp(fileType.c_str(), "CC  ") == 0) {
		cout<<"Loading Binary Case Control Data\n";
		theSample = new BasicSample(0, 0, 0, 0.0, 0.0, missingData, "frombin");
	} 
	else {
		cout<<"Unrecognized file type: "<<fileType<<". Unable to continue\n";
		return NULL;
	}

	if (theSample) {
		theSample->LoadBinarySample(genotypes, metadata, peopleCount, chromCount, chromSizes);
	}
	
	return theSample;	
}

void Sample::WriteBinaryHeader(ostream &genotypes, uint fileVersion, const char *type){
	//Header for binary format:
	//   SIZE           	//    Description       //     VALUE														//
	//   4 char  			//    File Type         //     PED															//
	//   2 char          //    bitsex            //     LE, BE (little endian, big endian)					//
 	//	  4 byte float   	// 	version				//     1.0															//
   //   4 byte uint    	//    bits per block    //     (for now, whatever dynamicbitset defaults to	//
	//   4 byte uint    	//    flags					//     For now, 1 means there is missing data, 0 means none
	//   4 byte uint    	//    chromosome count  //
   //   4 byes / chrom  //    loci count			//
 	
	uint chromCount = people[0]->GetChromosomeCount();
	uint flags = 0x0;
	if (missingData > 0.0) 
		flags = flags|1;
	
   genotypes.write(type, 4);
	genotypes.write("LE", 2);
	genotypes.write((const char *)&fileVersion, 4);
	uint bpb = boost::dynamic_bitset<>::bits_per_block;
	genotypes.write((const char *)&bpb, 4);
	genotypes.write((const char *)&flags, 4);

	uint peopleCount = people.size();
	genotypes.write((const char *)&peopleCount, 4);
	genotypes.write((const char *)&chromCount, 4);

	for (uint i=0; i<chromCount; i++) {
		uint c = people[0]->GetSnpCount(i);
		genotypes.write((const char *)&c, 4);
	}

}
/*
int BasicSample::LoadBinaryDataset(ostream &metadata, ostream &genotypes) {
	uint count = LoadBinaryHeader(genotypes);

	people.clear();

	for (uint i=0; i<count; i++) {
		Individual person;
		person.LoadMDR(metadata, genotypes);
		people.push_back(person);
	}
	return count;
}
*/
int BasicSample::WriteBinaryDataset(ostream &metadata, ostream &genotypes) {
	WriteBinaryHeader(genotypes, 1, "CC  ");
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WriteMDR(metadata, genotypes, missingData > 0.0, affCount > 0);
	return count;
}

void BasicSample::WriteSnpDetails(const char *project, PoolManager *pools) {
	stringstream filename;
	filename<<project<<".map";
	ofstream mapFile(filename.str().c_str());
	PoolManager::Iterator itr = pools->GetIterator();
	ChromPool *pool = itr.GetNext();

	while (pool) {
		LocusArray &loci = pool->GetLoci();
		LocusArray::iterator loc = loci.begin();
		LocusArray::iterator end = loci.end();

		while (loc != end) {
			mapFile<<pool->GetLabel()<<" "<<loc->GetLabel()<<" "<<loc->GetLocation()<<"\n";
			loc++;
		}
		pool = itr.GetNext();
	}
}
void BasicSample::LoadBinarySample(ifstream *genotypes, ifstream *meta, uint peopleCount, uint chromCount, vector<uint> &chrID) {
	for (uint i=0; i<peopleCount; i++) {
		Individual *ind = new Individual(i, 1000, chromCount);
		ind->ReadBinaryMdr(genotypes, missingData > 0.0, true, chrID);
		people.push_back(ind);
	}
}

bool ContinuousSample::Verify(PenetranceModel* model) {
	bool isValid = model->IsContinuous();
	if (isValid) {
		float lowerTailMax, upperTailMin;
		((ContinuousModel*)model)->GetTailBounds(lowerTailMax, upperTailMin);
		isValid = lowerTailMax + upperTailMin + midCount + unaffCount == 0 || (lowerTailMax + upperTailMin != 0 && midCount + unaffCount > 0);
		if (!isValid) {
			cerr<<"Unable to create a continuous model with tails without specifying a proper grand mean. Please see manual for instructions on setting the grand mean for the disease model.\n";
		}
	}
	else {
		cerr<<"Users must define a continuous model to be used with continuous datasets. Please see the manual for instructions on properly configuring continuous datasets. \n";
	}
	return isValid;
}
/**
 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
 * @param os The stream to which the sample will be written	

void BasicSample::WritePhased(ostream &os) {
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WritePhased(os, i+1);

} */
void BasicSample::GenerateReport(ostream &os, uint padding) {
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(genoError*100.0)<<"% Genotype Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(phenoError*100.0)<<"% Phenocopy Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(missingData*100.0)<<"% Missing Data "<<endl;
	os<<setw(padding - 5)<<"Case Control Sample : "<<affCount<<"A/"<<unaffCount<<"U"<<endl;
}


void ContinuousSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model) {
	uint sampleSize = people.size();
	uint errorCount = (uint)((float)affCount * phenoError);
	vector<int> indices;
	for (uint i=0; i<sampleSize; i++) 
		indices.push_back(i);
	random_shuffle(indices.begin(), indices.end(), Utility::Random::globalGenerator);

	BitSetType tested(sampleSize, false);
	Individual *individual=NULL;
	int n=0;
	for (uint i=0; i<errorCount; i++) {
		int idx = indices[n++];
		if (!tested[idx]) {
			individual = people[idx];
			float outcome;
			int status = 3;
			
			int id = individual->GetID();
			do {
				delete individual;
				individual = new Individual(id, id, pools->GetPoolCount());
				PoolManager::Iterator itr = pools->GetIterator();
				ChromPool *pool = itr.GetNext();
	
				while (pool) {
					pool->ApplyChromosomeIDs(individual);
					//At this point, the pools may or may not be in memory
					pool->ApplyPresentGenotypes(individual, true);
					pool=itr.GetNext();
				}				
				status = ((ContinuousModel*)model)->DrawPopulationStatus(outcome);
			} while (status>2);
			
			individual->SetStatus(status);
			people[idx]=individual;
			//Need to make sure that the tails are handled properly. In otherwords, if an individual
			//falls into a region of inviable outcome, we need to delete them and start over again
			individual->SetOutcome(outcome);
			tested[idx]=true;
		}
	}
}


}
