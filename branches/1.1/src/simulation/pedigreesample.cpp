//
// C++ Implementation: pedigreesample
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreesample.h"
#include "poolmanager.h"
#include <iomanip>
#include "meantable.h"

#define FILEVERSION 1.0
#define HAS_MISSING 


namespace Simulation {

bool PedigreeSample::UseOriginalCross = false;

using namespace Utility;

PedigreeSample::PedigreeSample(float genoError, float phenoError, float missing, const char *desc) : 
	Sample(genoError, phenoError, missing), pedID(0), familyCount(0), description(desc) { }

void PedigreeSample::AddFamilyType(uint aff, uint unaff, uint count, uint extras) {
	FamilyMakeup newFam(aff, unaff, count, extras);
	familyTypes.push_back(newFam);
	familyCount += count;
}

PedigreeSample::~PedigreeSample() { }

void PedigreeSample::ResetPedID() { 
	pedID=1;
}


void PedigreeSample::ApplyPresentGenotypes(PoolManager &pools, Individual &person) {
	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();

	while (pool) {
		pool->ApplyPresentGenotypes(&person, true);
		pool = itr.GetNext();
	}
		
}

void PedigreeSample::ReserveIndividual(PoolManager &pools, Individual &person) {
	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();
	while (pool) {
		pool->ReserveIndividual( &person);
		pool  = itr.GetNext();
	}
}	


//With SIMLA, we no longer need other forms of heterogeneity
//void PedigreeSample::BuildSample(PoolManager &pools, ModelManager& models, uint familyCount) {
void PedigreeSample::BuildSample(PoolManager &pools, PenetranceModel *model) {
	affectedChildren =0;
	Purge();
	familyStartPositions.clear();

//	PenetranceModel *model;
	uint count = familyTypes.size();
	for (uint i = 0; i<count; i++) {
		FamilyMakeup &fam = familyTypes[i];
		BuildSample(pools, model, fam);
	}
}

void PedigreeSample::BuildSample(PoolManager &pools, PenetranceModel* model, FamilyMakeup& family) {
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models	
	uint minChildCount = family.AffectedSibs() + family.UnaffectedSibs();

	uint affectedSibCount = family.AffectedSibs();		///<Used in determining what a valid pedigree is
	uint unaffectedSibCount = family.UnaffectedSibs();	///<Used to determine what a valid pedigree is
	int extraCount = family.ExtraSibs();
	uint individualCount = family.GetFamilyCount();
	family.ResetPedigrees();
	while (individualCount-- > 0) {
		pedID++;
		bool validPedigree=false;

		uint attempts = 0;
		
		while (!validPedigree) {
			uint childCount = minChildCount + Utility::Random::globalGenerator(extraCount + 1);

			if (attempts++ > 10000000) {
				cout<<"Unable to find a complete dataset after 10000000 tries\n";
				attempts=0;
			}
			uint id=1;

			uint affecteds = 0;
			Individual *dad = DrawIndividual(id++, pedID, &pools, model, -1, false);
			Individual *mom = DrawIndividual(id++, pedID, &pools, model, -1, true);
			

			Individual **children = new Individual*[childCount];

			//Build up the children
			for (uint i=0; i<childCount; i++) {
				Individual *child = children[i];
					child=mom->DLCross(dad, id++);
				children[i]= child;
				ApplyPresentGenotypes(pools, *child);
				if (child->ApplyStatus(model))  {
//					cout<<"\nX";
					affecteds++;
				} 
			}

			

			//Test for appropriate familial shape
			if (affecteds >= affectedSibCount && childCount - affecteds >= unaffectedSibCount) {
				affectedChildren+=affectedSibCount;
				family.AddPedigree(pedID);
				validPedigree=true;
				//Record the start position in our special array
				FamilyPositionNode pNode(people.size(), childCount+2, affecteds, childCount - affecteds);
				familyStartPositions.push_back(pNode);

				people.push_back(dad);
				people.push_back(mom);
				affecteds=affectedSibCount;
				for (uint i=0; i<childCount; i++) {	
					people.push_back(children[i]);
					children[i] = NULL;
				}
				ReserveIndividual(pools, *dad);
				ReserveIndividual(pools, *mom);
			}
			else {
				delete dad;
				delete mom;
			}
			//OK, let's clean up any children that didn't go into datasets, if there are any
			for (uint i=0; i<childCount; i++) 
				if (children[i])
					delete children[i];
			delete[] children;
		}
	}

}

/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 */
int PedigreeSample::WriteDataset(ostream &os, uint *gtCount) {

	vector<FamilyPositionNode> startPoints = familyStartPositions;
	random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	
	uint count = startPoints.size();

	uint pos = 0;
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		pos = startPoints[i].start;

		uint familySize = startPoints[i].count;

		for (uint o=0;o<familySize; o++) {
			Individual *person=people[pos++];
			if (person->DoIncludeInDataset()) {
				observations++;
				person->WritePedigree(os, gtCount);
			}
		}
	}
	return observations;
}



int PedigreeSample::WriteBinaryDataset(ostream &meta, ostream &genotypes) {
	uint observations = 0;
	uint typeCount = familyTypes.size();
	for (uint i = 0; i<typeCount; i++) 
		observations += WriteBinaryDataset(meta, genotypes, familyTypes[i]);
	
	return observations;
}

void PedigreeSample::WriteSnpDetails(const char *project, PoolManager *pools) {
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
			mapFile<<pool->GetLabel()<<" "<<loc->GetLabel()<<" "<<loc->MapPosition()<<"\n";
			loc++;
		}
		pool = itr.GetNext();
	}
	mapFile.close();

	//Now the DAT file, for merlin users
	filename.str("");
	filename<<project<<".dat";
	mapFile.open(filename.str().c_str());
	itr = pools->GetIterator();
	pool = itr.GetNext();
	mapFile<<"A Disease\n";

	while (pool) {
		LocusArray &loci = pool->GetLoci();
		LocusArray::iterator loc = loci.begin();
		LocusArray::iterator end = loci.end();

		while (loc != end) {
			mapFile<<"M "<<loc->GetLabel()<<"\n";
			loc++;
		}
		pool = itr.GetNext();
	}
	mapFile.close();
}





void PedigreeSample::LoadBinarySample(ifstream *genotypes, ifstream *meta, uint peopleCount, uint chromCount, vector<uint> &chrID) {
	cout<<"We haven't done binary pedigree loading yet!\n";
	assert(0);
}


int PedigreeSample::WriteDataset(ostream &os, uint *gtCount, FamilyMakeup& fam) {
	
	uint count = people.size();


	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		Individual *person=people[i];
		if (person->DoIncludeInDataset()) {
			observations++;
			person->WritePedigree(os, gtCount);
		}
	}
	return observations;

}

int PedigreeSample::WriteBinaryDataset(ostream &meta, ostream &genotypes, FamilyMakeup& fam) {
	WriteBinaryHeader(genotypes, 1, "PED ");

	vector<FamilyPositionNode> startPoints = familyStartPositions;
	random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	uint count = startPoints.size();

	uint pos = 0;
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		pos = startPoints[i].start;
		uint familySize = startPoints[i].count;

		for (uint o=0;o<familySize; o++) {
			Individual *person=people[pos++];
			if (person->DoIncludeInDataset()) {
				observations++;
				person->WritePedigree(meta, genotypes, missingData > 0.0);
			}
		}
	}
	return observations;
}


/** I need to figure out how to id the mom and dad 
void PedigreeSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model) {
	uint famTypeCount = familyTypes.size();
	for (uint i=0; i<famTypeCount; i++) 
		ApplyPhenocopyError(pools, familyTypes[i], model);
}
*/

vector<Individual*> *PedigreeSample::BuildPedigree(PoolManager *pools, PenetranceModel* model, uint pedID, uint affected, uint unaffected, int extra) {
	int id = 1;
	int childCount = affected+unaffected+extra;
	vector<Individual *> *individuals = new vector<Individual*>;
	bool doContinue = true;
	uint affecteds = 0;
	while (doContinue) {
		Individual*dad = DrawIndividual(id++, pedID, pools, model, -1, false);
		Individual *mom = DrawIndividual(id++, pedID, pools, model, -1, true);
		individuals->push_back(dad);
		individuals->push_back(mom);
		affecteds = 0;
		
		for (int i=0; i<childCount; i++) {
			Individual *child = dad->DLCross(mom, id++);
			individuals->push_back(child);
			ApplyPresentGenotypes(*pools, *child);
			if (child->ApplyStatus(model) == 1)
				affecteds++;
		}
		doContinue = !(affecteds >= affected && childCount - affecteds >= unaffected);
	
		if (doContinue) {
			for (int i=0; i<2+childCount; i++) 
				delete (*individuals)[i];
			individuals->clear();
		}

	}
	return individuals;
}


void PedigreeSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model) {
	//This structure has the location and size of each pedigree
	vector<FamilyPositionNode*> pedigrees;
	vector<FamilyPositionNode>::iterator familyStart = familyStartPositions.begin();
	vector<FamilyPositionNode>::iterator familyEnd = familyStartPositions.end();
	while (familyStart != familyEnd) {
		familyStart->phenocopies=0;
		FamilyPositionNode *fp = &(*familyStart);
		for (uint i=0; i<familyStart->affecteds; i++) {
			pedigrees.push_back(fp);
		}
		familyStart++;
	}

	random_shuffle(pedigrees.begin(), pedigrees.end(), Utility::Random::globalGenerator);
	
	int phenocopies = (int)((float)affectedChildren*phenoError);
	for (int i=0; i<phenocopies; i++) 
		pedigrees[i]->phenocopies++;
		
	
	
	vector<FamilyPositionNode*>::iterator itr = pedigrees.begin();
	vector<FamilyPositionNode*>::iterator end = pedigrees.end();

	while (itr != end) {

		FamilyPositionNode *node = *itr;
		int members = node->count;
		//For now, let's push through until we run out of phenocopies
		int start = node->start;
		int phenoCount = node->phenocopies;
		if (phenoCount>0) {
			vector<Individual*>* newPed = BuildPedigree(pools, model, people[start]->GetPedigreeID(), node->affecteds - phenoCount, node->unaffecteds + phenoCount, 0);
			delete people[start];
			delete people[start+1];
			people[start] = (*newPed)[0];
			people[start+1] = (*newPed)[1];

			for (int i=2; i<members; i++) {
				delete people[start+i];
				Individual *ind = newPed->at(i);
				people[start+i]=ind;
				if (phenoCount > 0 && !ind->IsAffected()) {
					phenoCount--;
					ind->SetStatus(1);
				}
			}
		}
		itr++;
	}

}



void PedigreeSample::GenerateReport(ostream &os, uint padding) {
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(genoError*100.0)<<"% Genotype Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(phenoError*100.0)<<"% Phenocopy Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(missingData*100.0)<<"% Missing Data "<<endl;

	uint famTypeCount = familyTypes.size();
	for (uint i=0; i<famTypeCount; i++) 
		os<<setw(padding)<<"Pedigree Sample "<<i + 1<<"("
		  <<familyTypes[i].GetFamilyCount()<<") : "	<<familyTypes[i].AffectedSibs()
		  <<"A/"<<familyTypes[i].UnaffectedSibs()<<"U "<<endl;

}



void PedigreeMixedSample::BuildSample(PoolManager &pools, PenetranceModel* model) {
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models
	abort();

	Purge();
	uint individualCount = familyCount;
	familyStartPositions.clear();
	while (individualCount-- > 0) {
		pedID++;
		bool validPedigree=false;

		uint attempts = 0;

		while (!validPedigree) {

			//Random number generator returns a value from 0 to n-1
			uint childCount = Utility::Random::globalGenerator((int)maxChildren) + 1;

			if (attempts++ > 10000) {
				cout<<"Unable to find a complete dataset after 10000 tries\n";
				attempts=0;
			}
			uint id=1;

			uint affecteds = 0;


			Individual *dad = DrawIndividual(id++, pedID, &pools, model, -1, false);
			Individual *mom = DrawIndividual(id++, pedID, &pools, model, -1, true);


			Individual **children = new Individual*[maxChildren];

/*			Individual *dad = new Individual(id++, pedID, chromCount);
			Individual *mom = new Individual(id++, pedID, chromCount);
			DrawIndividual(pools, *mom);
			DrawIndividual(pools, *dad);
*/
			//Build up the children
			for (uint i=0; i<maxChildren; i++) {
				//children[i]->Init(id++, pedID, chromCount);
				//This is broken! I am not using crossover to generate children!
				children[i]=dad->Cross(mom, id++);
				
				if (i < childCount) {
					if (children[i]->ApplyStatus(model))  {
						affecteds++;
					}
				}
				else 		//This will let us skip this child
					children[i]->DoIncludeInDataset(false);
			}

			//Test for appropriate familial shape
			if (affecteds > 0) {
				validPedigree=true;
				//Record the start position in our special array
				familyStartPositions.push_back(FamilyPositionNode(people.size(), maxChildren, affecteds, maxChildren-affecteds) );
				people.push_back(dad);
				people.push_back(mom);
				ReserveIndividual(pools, *dad);
				ReserveIndividual(pools, *mom);
				for (uint i=0; i<maxChildren; i++) {	
					people.push_back(children[i]);
				}
			}
			else {
				delete dad;
				delete mom;
				for (uint i=0; i<maxChildren; i++) 
					delete children[i];

			}
			delete[] children;
		}
	}

}

string PedigreeSample::GetSummary() { 
	vector<FamilyMakeup>::iterator itr = familyTypes.begin();
	vector<FamilyMakeup>::iterator end = familyTypes.end();

	int count = 0;
	while (itr != end) {
		int contribution = (2+itr->AffectedSibs()+itr->UnaffectedSibs()+(itr->ExtraSibs()/2));
		count+=(contribution*itr->GetFamilyCount());
		itr++;
	}
	stringstream ss;
	ss<<familyTypes.size()<<" Family Types\t~"<<count<<" Individuals";
	return ss.str();

}

string PedigreeSample::GetDetails() { 
	vector<FamilyMakeup>::iterator itr = familyTypes.begin();
	vector<FamilyMakeup>::iterator end = familyTypes.end();
	stringstream ss;
	ss<<"DATASET PED "<<description<<" "<<genoError<<" "<<phenoError<<" "<<missingData<<"\n";
	
	while (itr != end) {
		ss<<itr->GetConfigurationDetails()<<"\n";
		itr++;
	}
	return ss.str();

}

void PedigreeContinuousSample::WriteDatFile(PoolManager *pools, ostream& os) {
	os<<"A Disease\n";
	os<<"T outcome\n";
	PoolManager::Iterator itr = pools->GetIterator();
	ChromPool *pool = itr.GetNext();

	while (pool) {
		LocusArray &loci = pool->GetLoci();
		LocusArray::iterator loc = loci.begin();
		LocusArray::iterator end = loci.end();

		while (loc != end) {
			os<<"M "<<loc->GetID()<<"\n";
			loc++;
		}
	}

}


bool PedigreeContinuousSample::Verify(PenetranceModel* model) {
	bool isValid = model->IsContinuous();
	if (isValid) {
/*		float lowerTailMax, upperTailMin;
		((ContinuousModel*)model)->GetTailBounds(lowerTailMax, upperTailMin);
		isValid = lowerTailMax + upperTailMin + midCount + unaffCount == 0 || (lowerTailMax + upperTailMin != 0 && midCount + unaffCount > 0);
		if (!isValid) {
			cerr<<"Unable to create a continuous model with tails without specifying a proper grand mean. Please see manual for instructions on setting the grand mean for the disease model.\n";
		}
 */
	}
	else {
		cerr<<"Users must define a continuous model to be used with continuous datasets. Please see the manual for instructions on properly configuring continuous datasets. \n";
	}
	return isValid;
}

/** I need to figure out how to id the mom and dad */
void PedigreeContinuousSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model) {
	uint famTypeCount = familyTypes.size();
	for (uint i=0; i<famTypeCount; i++) 
		ApplyPhenocopyError(pools, familyTypes[i], model);
}
void PedigreeContinuousSample::ApplyPhenocopyError(PoolManager *pools, FamilyMakeup& fam, PenetranceModel* model) {
	uint sampleSize = people.size();
	uint affectedSibCount = fam.AffectedSibs() + fam.UnaffectedSibs() + fam.ExtraSibs();
	uint errorCount = (uint)((float)(familyCount * affectedSibCount) * phenoError);

	BitSetType tested(sampleSize, false);

	Individual *individual=NULL;
	for (uint i=0; i<errorCount; i++) {
		bool canBeUsed = false;
		uint idx;
		while (!canBeUsed) {
			idx = Utility::Random::globalGenerator((int)sampleSize);
			if (!tested[idx]) {
				individual = people[idx];
				canBeUsed = individual->GetMother() && individual->GetFather();
				tested[idx]=true;
			}
		}
		Individual *mom = individual->GetMother();
		Individual *dad = individual->GetFather();
			
		uint newID = individual->GetID();
		int status = 3;
		float outcome;

		Individual *newIndividual = NULL;
		do {
			newIndividual = dad->Cross(mom, newID);
			status = ((StatusModel::ContinuousModel*)model)->DrawPopulationStatus(outcome);
			//If it's inviable, we have to start over.
			if (status == 3) {
				delete newIndividual;
			}
		} while (status == 3);
		delete individual;
		newIndividual->SetOutcome(outcome);
		newIndividual->SetStatus( status );
		people[idx]=newIndividual;
	}
}


int PedigreeContinuousSample::WriteDataset(ostream &os, uint *gtCount, FamilyMakeup& fam) {
	
	uint count = people.size();


	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		Individual *person=people[i];
		if (person->DoIncludeInDataset()) {
			observations++;
			person->WritePedigree(os, gtCount);
		}
	}
	return observations;

}
int PedigreeContinuousSample::WriteBinaryDataset(ostream &meta, ostream &genotypes) {
	uint observations = 0;
	uint typeCount = familyTypes.size();
	for (uint i = 0; i<typeCount; i++) 
		observations += WriteBinaryDataset(meta, genotypes, familyTypes[i]);
	
	return observations;
}
int PedigreeContinuousSample::WriteBinaryDataset(ostream &meta, ostream &genotypes, FamilyMakeup& fam) {
	
	WriteBinaryHeader(genotypes, 1, "PED ");

	vector<FamilyPositionNode> startPoints = familyStartPositions;
	random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	
	uint count = startPoints.size();


	uint pos = 0;
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		pos = startPoints[i].start;

		uint familySize = startPoints[i].count;

		for (uint o=0;o<familySize; o++) {
			Individual *person=people[pos++];
			if (person->DoIncludeInDataset()) {
				observations++;
				person->WritePedigree(meta, genotypes, missingData > 0.0);
			}
		}
	}
	return observations;

	return 0;
}
/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 */
int PedigreeContinuousSample::WriteDataset(ostream &os, uint *gtCount) {

	vector<FamilyPositionNode> startPoints = familyStartPositions;
	random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	
	uint count = startPoints.size();

	uint pos = 0;
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		pos = startPoints[i].start;

		uint familySize = startPoints[i].count;

		for (uint o=0;o<familySize; o++) {
			Individual *person=people[pos++];
			if (person->DoIncludeInDataset()) {
				observations++;
				person->WritePedigree(os, gtCount, true);
			}
		}
	}
	return observations;
}
void PedigreeContinuousSample::AppendGenotypeCountsToReport(StatusModel::ModelLociArray loci) {
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
int PedigreeSample::WriteDataset(const char *filename, uint *gtCount) {
	ofstream file(filename);
	uint observations = WriteDataset(file, gtCount);
	file.close();	
	return observations;
}

void PedigreeContinuousSample::ReportGenotypeCounts(vector<StatusModel::DiseaseLocus> loci, ostream& output) {
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


void PedigreeContinuousSample::WriteSnpDetails(const char *project, PoolManager *pools) {
	stringstream filename;
	filename<<project<<".dat";
	ofstream datFile(filename.str().c_str());
	PoolManager::Iterator itr = pools->GetIterator();
	ChromPool *pool = itr.GetNext();
	datFile<<"A Disease\n";
	datFile<<"T Outcome\n";

	while (pool) {
		LocusArray &loci = pool->GetLoci();
		LocusArray::iterator loc = loci.begin();
		LocusArray::iterator end = loci.end();

		while (loc != end) {
			datFile<<"M "<<loc->GetLabel()<<"\n";
			loc++;
		}
		pool = itr.GetNext();
	}
}

void PedigreeContinuousSample::BuildSample(
			PoolManager &pools, PenetranceModel* model, FamilyMakeup& family) {
	//OK, determine the basic size of the family
	uint minChildCount = family.AffectedSibs() + family.UnaffectedSibs() + family.ExtraSibs();
	uint lowerTail 	= family.AffectedSibs();		///<Used in determining what a valid pedigree is
	uint upperTail 	= family.UnaffectedSibs();		///<Used to determine what a valid pedigree is
	uint mids 		= family.ExtraSibs();

	float lowerThresh, upperThresh;
	((ContinuousModel*)model)->GetTailBounds(lowerThresh, upperThresh);
	if (lowerTail == upperTail) {
		lowerTail = minChildCount;
		upperTail = mids = 0;
	}

	uint individualCount = family.GetFamilyCount();

	while (individualCount-- > 0) {
		pedID++;
		bool validPedigree=false;

		uint attempts = 0;
		
		while (!validPedigree) {
			uint childCount = minChildCount;// + Utility::Random::globalGenerator(extraCount + 1);

			if (attempts++ > 10000000) {
				cout<<"Unable to find a complete dataset after 10000000 tries\n";
				attempts=0;
			}
			uint id=1;

			Individual *dad = DrawIndividual(id++, pedID, &pools, model, -1, false);
			Individual *mom = DrawIndividual(id++, pedID, &pools, model, -1, true);

			Individual **children = new Individual*[childCount];

			vector<uint> stati(3, 0);

			//Build up the children. For the time being, we call any viable offspring as affected
			for (uint i=0; i<childCount; i++) {
				Individual *child = children[i];
				child=dad->DLCross(mom, id++);
				children[i]= child;
				ApplyPresentGenotypes(pools, *child);
				int status = child->ApplyStatus(model);
				if (status == 3)
					i--;
				else 
					stati[status]++;
			}

			if (stati[0] == lowerTail && stati[1] == upperTail && stati[2] == mids) {
				//Test for appropriate familial shape	
				validPedigree=true;
				//Record the start position in our special array
				FamilyPositionNode pNode(people.size(), childCount+2, stati[0]+stati[1], stati[2]);
				familyStartPositions.push_back(pNode);
	
				people.push_back(dad);
				people.push_back(mom);
				for (uint i=0; i<childCount; i++) {	
					people.push_back(children[i]);
					children[i] = NULL;
				}
				ReserveIndividual(pools, *dad);
				ReserveIndividual(pools, *mom);
			}
			else {
				delete mom;
				delete dad;
			}
			delete[] children;
		}
	}

}


}
