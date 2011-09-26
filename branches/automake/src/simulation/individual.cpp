//
// C++ Implementation: individual
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sstream>
#include "individual.h"
#include "penetrancemodel.h"
namespace Simulation {

using namespace std;
using namespace Utility;
using namespace StatusModel;

/**
 * @brief Standard pedigree is treated as the 10 column version
 */
bool Individual::StandardPedigreeHeader 	= false;
bool Individual::PhasedPedigrees			= false;
float Individual::FtoM_BirthRatio			= 0.47619;

/**
 * @brief If an allele's minor allele frequency is lower than this, we won't write it to datasets 
 */
float Individual::minAlFreqThreshold		= 0.00;

Individual::Individual(uint id, uint pedID, uint chrCount /*=1*/) : id(id), pedID(pedID), patID(0), 
			matID(0), includeInDataset(true), status(0), gender(0), chromosomeCount(chrCount), mom(NULL), 
			dad(NULL), chrID1(chrCount), chrID2(chrCount)
#ifdef USE_XY
			,par(NULL)
#endif

{
	matChrom = new Chromosome[chrCount];
	patChrom = new Chromosome[chrCount];
	genePoolIDs = new string[chrCount];
	missingData = new BitSetType[chrCount];
}

Individual::Individual() : id(0), pedID(0), patID(0), matID(0), includeInDataset(true), status(0), gender(0),
		 matChrom(NULL), patChrom(NULL), chromosomeCount(0), genePoolIDs(NULL), mom(NULL), 
		dad(NULL), chrID1(0), chrID2(0) 
#ifdef USE_XY
			,par(NULL)
#endif

{ 
}

void Individual::Init(uint id, uint pedID, uint chrCount /*=1*/)  {
	id = id;
	pedID = pedID;
	chromosomeCount = chrCount;
	if (matChrom) 
		matChrom = new Chromosome[chrCount];

	if (patChrom)
		patChrom = new Chromosome[chrCount];
	if (genePoolIDs)
		delete[] genePoolIDs;
	genePoolIDs = new string[chrCount];
	chrID1.resize(chrCount);
	chrID2.resize(chrCount);
}

bool Individual::DoIncludeInDataset() {
	return includeInDataset;
}

void Individual::DoIncludeInDataset( bool doinclude){
	includeInDataset=doinclude;
}

float Individual::MinorAlleleFreq(uint chID, uint locID) {
	float maf = 0.0;				///<Minor allele freq
	if (patChrom) {
		Locus l = patChrom[chID].GetLocus(locID);
		maf = l.GetMinAlleleFreq();
	}
	return maf;		
}

Individual::~Individual()	{
	if (matChrom)
		delete[] matChrom;
	if (patChrom)
		delete[] patChrom;
	if (genePoolIDs)
		delete[] genePoolIDs;
	if (missingData)
		delete[] missingData;
}




void Individual::SetDLPedigreeMeta(Individual *dad, Individual *mom) {
	this->patID=dad->GetID();
	this->matID=mom->GetID();
	this->dad=dad;
	this->mom=mom;
#ifdef USE_XY
	if (dad->chromXY.chrom1) {
		bool isFemale = Random::globalGenerator.drand() < FtoM_BirthRatio;
		if (isFemale) {
			AlleleSource<LocusXY> *d = dad->chromXY.Draw(Random::globalGenerator, dad->par, isFemale);
			assert(VerifyChromosomeIs(d, true));
			assert(dad->VerifyChromosomeIs(false));
			AlleleSource<LocusXY> *m = mom->chromXY.Draw(Random::globalGenerator, NULL, isFemale);
			assert(mom->VerifyChromosomeIs(true));
			assert(VerifyChromosomeIs(m, true));
			SetXX(m, d, dad->par);
		}
		else {
			AlleleSource<LocusXY> *d = dad->chromXY.Draw(Random::globalGenerator, dad->par, isFemale);
			assert(VerifyChromosomeIs(d, false));
			assert(dad->VerifyChromosomeIs(false));
			AlleleSource<LocusXY> *m = mom->chromXY.Draw(Random::globalGenerator, NULL, isFemale);
			assert(mom->VerifyChromosomeIs(true));
			assert(VerifyChromosomeIs(m, true));
			SetXY(m, d, dad->par);
		}
	}
#endif
	SetChromosomalData(dad->DLCross(), mom->DLCross());
}

void Individual::SetPedigreeMeta(Individual *dad, Individual *mom) {
	this->patID=dad->GetID();
	this->matID=mom->GetID();
	this->dad=dad;
	this->mom=mom;
	SetChromosomalData(dad->Cross(), mom->Cross());
}

void Individual::SetPoolID(uint chromeID, const char *id) {
	assert(chromeID < chromosomeCount && genePoolIDs);
	genePoolIDs[chromeID]=id;
}

std::string Individual::GetPoolID(uint chromID) {
	assert(chromID <chromosomeCount);
	return genePoolIDs[chromID];
}

void Individual::SetChromosomalData(Chromosome *p, Chromosome *m) {
	if (matChrom)
		delete[] matChrom;
	if (patChrom)
		delete[] patChrom;

	matChrom = m;
	patChrom = p;
	
	int chromCount = chromosomeCount;
#ifdef USE_XY
	if (chromXY.chrom1) 
		chromCount--;
#endif
	for (int i=0; i<chromCount; i++)
		missingData[i].resize(p[i].LociCount(), false);
}

void Individual::SetChromIDs(uint chromID, uint chID1, uint chID2) {
	chrID1[chromID] = chID1; 
	chrID2[chromID] = chID2;
}

bool Individual::NeedsCompleteGenotypeInformation() {
	return mom == NULL;
}

void Individual::SetChromosomalData(uint chromID, Chromosome &p, Chromosome &m) {
	assert(chromID < chromosomeCount);
	matChrom[chromID] = m;
	patChrom[chromID] = p;
	missingData[chromID].resize(p.LociCount(), false);
	//cout<<"Setting Missing ("<<chromID<<") to "<<p.LociCount()<<"\n";
}


uint Individual::GetStatus() { 
	return status;
}

void Individual::SetStatus(int isAffected) {
	status = isAffected;
}

void Individual::SetGender(int gender) {
	this->gender = gender;
}

bool Individual::IsAffected() {
	return status != 0;
}


void Individual::SetOutcome(float& outcome) {
	this->outcome=outcome;
}
float Individual::GetOutcome() {
	return outcome;
}

bool Individual::MissingGenotype(uint chrID, uint locus) {
	return missingData[chrID][locus];
}
bool Individual::ClearLocus( uint chrID, uint locus) {
	assert(chromosomeCount > chrID);
	if (missingData[chrID].size() <= locus)
		cout<<"Hm, we are about to die! Missing data["<<chrID<<"].size() == "<<missingData[chrID].size()<<" & Locus == "<<locus<<"\n";
	assert(missingData[chrID].size() > locus);

	if (missingData[chrID][locus]) 
		return false;
	else 
		missingData[chrID][locus] = true;
	return true;
}

void Individual::InitAlleleSource() {
	string mchID = ToString((int)GetPedigreeID()) + ":" + ToString((int)id) + "m";
	string pchID = ToString((int)GetPedigreeID()) + ":" + ToString((int)id) + "p";
	for (uint i=0; i<chromosomeCount; i++) {
		matChrom[i].InitAlleleSource(mchID.c_str());
		patChrom[i].InitAlleleSource(pchID.c_str());
	}
}

Individual *Individual::Clone() {
	Individual *ind = new Individual(id, pedID, chromosomeCount);
	ind->patID = patID;
	ind->matID = matID;
	ind->includeInDataset = includeInDataset;
	ind->status = status;
	ind->outcome = outcome;
	if (matChrom) {
		//ind->matChrom = new Chromosome[chromosomeCount];
		for (uint i=0; i<chromosomeCount; i++)
			ind->matChrom[i] = matChrom[i];
	}
	if (patChrom) {
		//ind->patChrom = new Chromosome[chromosomeCount];
		for (uint i=0; i<chromosomeCount; i++) 
			ind->patChrom[i] = patChrom[i];
	}
#ifdef USE_XY
	ind->chromXY = chromXY;
#endif
	for (uint i=0; i<chromosomeCount; i++) {
		ind->genePoolIDs[i]=genePoolIDs[i];
		ind->missingData[i]=missingData[i];
	}
	if (mom)
		ind->mom = mom->Clone();
	if (dad)
		ind->dad = dad->Clone();

	ind->chrID1 = chrID1;
	ind->chrID2 = chrID2;
	ind->gender = gender;
	return ind;
}


/**
 * @Note errorDir == 0  -Allow for randomized direction
 */
bool Individual::ChangeGenotype(uint chrID, uint locus, int errorDir) {

	if (errorDir == 0)  {
		if (Utility::Random::globalGenerator((float)1.0) > 0.5)
			errorDir = 1;
		else
			errorDir = -1;
	}


	int gt = GetGenotype(chrID, locus) + errorDir;

	bool success=true;
	switch (gt) {
		case 0: 
			matChrom[chrID][locus] = false;
			patChrom[chrID][locus] = false;
			break;
		case 1:
			if (errorDir > 0) {
				if (Utility::Random::globalGenerator((float)1.0) > 0.5) 
					matChrom[chrID][locus]=true;
				else
					patChrom[chrID][locus]=true;
			}
			else {
				if (Utility::Random::globalGenerator((float)1.0) > 0.5) 
					matChrom[chrID][locus]=false;
				else
					patChrom[chrID][locus]=false;
			}
			break;
		case 2:
			matChrom[chrID][locus]=true;
			matChrom[chrID][locus]=true;
			break;
		default:
			success=false;
	}
	return success;
}

vector<uint> Individual::GetDiseaseGenotypes(DiseaseModel *model) {
	vector<uint> genotypes;
	if (model) {
		size_t modelSize = model->GetModelSize();
		for (uint i=0; i<modelSize; i++)  {
			StatusModel::DiseaseLocus l = model->GetLocus(i);
			genotypes.push_back(GetGenotype(l.chromosome, l.locusIdx));
		}
	}
	return genotypes;
}

int Individual::ApplyStatus(DiseaseModel *model) {
	int isAffected = 0;
	if (model) {
		vector<uint> genotypes = GetDiseaseGenotypes(model);

		//For continuous data, we need to know if there is any relevance here. If so, 
		//How do we determine the status relationship? And if it is irrelavent, how
		//do we show that it doesn't matter?
		isAffected = model->GetStatus(genotypes, outcome);
#ifdef DEBUG_VERBOSE
		cerr<<id<<":"<<pedID<<" [";
		for (int i=0; i<model->GetModelSize(); i++) 
			cerr<<" "<<genotypes[i];
		cerr<<" ] = "<<isAffected<<"\t";
		WriteMDR(cerr, NULL, true);
#endif //NDEBUG
	}
	else {
		if (Utility::Random::globalGenerator.drand() < 0.5) 
			isAffected = 1;
		outcome = 0.0;
	}
	SetStatus(isAffected);

	return isAffected;	
}



Individual *Individual::DLCross(Individual *father, size_t id, int gender/*= 0*/) {
	Individual *child = new Individual(id, father->GetPedigreeID(), chromosomeCount);
	child->gender = gender;
	child->SetDLPedigreeMeta(father, this);
	return child;
}

Individual *Individual::Cross(Individual *father, size_t id) {
//	Chromosome *pDNA = father->Cross();
//	Chromosome *mDNA = Cross();
	Individual *child = new Individual(id, father->GetPedigreeID(), chromosomeCount);
	child->SetPedigreeMeta(father, this);
//	child->SetChromosomalData(pDNA, mDNA);
	return child;
}

Chromosome *Individual::DLCross() {	
	int chromCount = chromosomeCount;
#ifdef USE_XY
	if (chromXY.chrom1)
		chromCount--;
#endif
	Chromosome *newChr = NULL;
	if (chromCount > 0) {
		newChr = new Chromosome[chromCount];
		size_t eventCount;
		for (int i=0; i<chromCount; i++)
			newChr[i] = patChrom[i].DLCross(matChrom[i], eventCount);
	}
	return newChr;
}

Chromosome *Individual::Cross() {
	Chromosome *newChr = new Chromosome[chromosomeCount];
	size_t eventCount;						///<need to capture the number of x-over events
	for (uint i=0; i<chromosomeCount; i++)
		newChr[i] = patChrom[i].Cross(matChrom[i], eventCount);

	return newChr;
}

int Individual::GetGenotype(uint chrID, uint locus) {

#ifdef USE_XY 
	if (chrID == (uint)-1)
		return chromXY.GetGenotype(locus);
#endif //USE_XY
	if (chromosomeCount <= chrID || locus >= matChrom[chrID].LociCount() || locus >= patChrom[chrID].LociCount()) {
		cout<<"-------------------------- "<<chrID<<" < "<<chromosomeCount<<" && "<<matChrom[chrID].LociCount()<<" < "<<locus<<" > "<<patChrom[chrID].LociCount()<<"\n";
		cout<<"--------------------------We asked for an invalid memory location! (individual.cpp:196)\n";
		assert(0);
	}
	return matChrom[chrID][locus] + patChrom[chrID][locus];
}
/*
uint Individual::GetChromosomeCount() {
	return chromosomeCount;
}

size_t Individual::GetSnpCount(uint chID) {
	return patChrom[chID].LociCount;
}*/
int Individual::GetSnpCount() {
	uint lociCount = 0;
	for (uint chID=0; chID<chromosomeCount; chID++) 
		lociCount+=patChrom[chID].LociCount();
#ifdef USE_XY
	if (chromXY.chrom1)
		lociCount+=chromXY.chrom1->GetLocusCount();
#endif //USE_XY
	return lociCount;
}

void Individual::ExpandModelLoci(uint chID, vector<uint>& modelLoci) {
	matChrom[chID].ExpandModelLoci(modelLoci);
	patChrom[chID].ExpandModelLoci(modelLoci);
}

void Individual::ResolveGenotypes(vector<uint> &modelLoci) {
	for (uint chID=0; chID<chromosomeCount; chID++) {
		ResolveGenotypes(chID, modelLoci);
	}
}


void Individual::PopulateChromosomePool(uint chID, std::vector<Chromosome>& pool) {
	pool.push_back(patChrom[chID]);
	pool.push_back(matChrom[chID]);
}

void Individual::ResolveGenotypes(uint chID, vector<uint> &modelLoci) {

		stringstream ss;
#ifdef DEBUG_VERBOSE
vector<uint> genotypes(modelLoci.size());
for (uint i=0; i<modelLoci.size(); i++) 
	genotypes[i] = GetGenotype(chID, modelLoci[i]);

WritePhased(cerr, id, false);
#endif
		patChrom[chID].ResolveGenotypes();
		matChrom[chID].ResolveGenotypes();	
#ifdef DEBUG_VERBOSE
WritePhased(cerr, id, false);

vector<uint> genotypesPost(modelLoci.size());
for (uint i=0; i<modelLoci.size(); i++)
	genotypesPost[i] = GetGenotype(chID, modelLoci[i]);

for (int i=0; i<modelLoci.size(); i++) {
	assert(genotypes[i] == genotypesPost[i]);
}
#endif
}

void Individual::ReferenceDistance(vector<double>& distances, int diseaseLocus) {
	//We don't know how to handle LOD for other chromosomes
	assert(chromosomeCount==1);
	matChrom[0].ReferenceDistance(distances, diseaseLocus);
}

void Individual::CalculateLOD(vector<int>& r, vector<int>& nr, int diseaseLocus) {
	//We don't want to calculate LOD for individuals that aren't born into the pedigree
	if (mom && dad) {
//WriteXOReport(cout);
		//We don't know how to handle LOD for other chromosomes
		assert(chromosomeCount==1);



		if (status == 1)  {
cerr<<"+";
//			matChrom[0].EvaluateRecombinants(r, diseaseLocus);
//			patChrom[0].EvaluateRecombinants(r, diseaseLocus);
//cerr<<"\n";
		}
		else  if (status == 0){
cerr<<"-";
//			matChrom[0].EvaluateRecombinants(nr, diseaseLocus);
//			patChrom[0].EvaluateRecombinants(nr, diseaseLocus);
//cerr<<"\n";
		}
	
		if (status < 2) {
			vector<int> recombinations(r.size());
			
			matChrom[0].EvaluateRecombinants(recombinations, diseaseLocus);
			patChrom[0].EvaluateRecombinants(recombinations, diseaseLocus);
	
			vector<int>::iterator itr = recombinations.begin();
			vector<int>::iterator end = recombinations.end();
			
			int idx =0;
			while (itr != end) {
				if (*itr++ > 0)
					r[idx++]++;
				else
					nr[idx++]++;
			}
		}
	}
//WritePhased(cerr, id, false);
}


void Individual::WriteForMendel(std::ostream& meta, std::ostream& genotypes) {
	static int statusConvertor[] = {1, 2, 0};
	meta<<pedID<<", ";								//Column #1
	meta<<id<<", ";									//Column #2
	meta<<matID<<", ";								//Column #3
	meta<<patID<<", ";								//Column #4
	meta<<gender<<", , ";								//Gender
	meta << statusConvertor[status] << "\n";


	boost::dynamic_bitset<unsigned char> data(GetSnpCount() * 2);
	int idx = 0;
	for (uint chID=0; chID<chromosomeCount; chID++)  {
		Chromosome &p = patChrom[chID];
		Chromosome &m = matChrom[chID];
		for (uint snp = 0; snp<p.LociCount(); snp++) {
			data[idx++]=p[snp];
			data[idx++]=m[snp];
		}
	}

#ifdef USE_XY
	cerr<<"Ooops, we haven't set mendel files up to work with XY chromosomes!\n";
	assert(0);
	if (chromXY.chrom1)
		chromXY.chrom1->WriteBinary(genotypes);
#endif //USE_XY
	uint totalBlocks = (data.size() / boost::dynamic_bitset<unsigned char>::bits_per_block) + 1;
assert(boost::dynamic_bitset<unsigned char>::bits_per_block == 8);
	dynamic_bitset<unsigned char>::block_type raw[totalBlocks];

	//Convert the bitset to raw data
	to_block_range(data, raw);

	//cout<<"writing "<<totalBlocks<<" ints to the stream "<<(sizeof(dynamic_bitset<>::block_type)*totalBlocks)<<"\n";
	genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	


}

void Individual::WritePedigree(std::ostream& meta, std::ostream&genotypes, bool writeMissingData, bool doWriteStatus /*=true*/) {
	static int statusConvertor[] = {1, 2, 0};
	meta<<pedID<<" ";								//Column #1
	meta<<id<<" ";								//Column #2
	meta<<matID<<" ";								//Column #3
	meta<<patID<<" ";								//Column #4

	if (StandardPedigreeHeader) {
		meta<<"0 ";			//uint nextOffspring;	//Column #5
		meta<<"0 ";			//uint nextPatSib;		//Column #6
		meta<<"0 ";			//uint nextMatSib;		//Column #7
	}
	if (id<3)
		meta<<id<<" ";
	else
		meta<<"0 ";			//uint sex;				//Column #8
	if (StandardPedigreeHeader) {
		meta<<"0 ";			//uint probStat;		//Column #9
	}
	meta << statusConvertor[status] << " ";

	for (uint chID=0; chID<chromosomeCount; chID++) 
		patChrom[chID].WriteBinary(genotypes);

#ifdef USE_XY
	if (chromXY.chrom1)
		chromXY.chrom1->WriteBinary(genotypes);
#endif //USE_XY

	for (uint chID=0; chID<chromosomeCount; chID++) 
		matChrom[chID].WriteBinary(genotypes);
	
	if (writeMissingData) {
		for (uint chID=0; chID<chromosomeCount; chID++) {
			BitSetType &misses = missingData[chID];
			uint totalBlocks = (misses.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
			dynamic_bitset<>::block_type raw[totalBlocks];

			//Convert the bitset to raw data
			to_block_range(misses, raw);

			genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	
		}
	}
#ifdef USE_XY
	if (chromXY.chrom2)
		chromXY.chrom2->WriteBinary(genotypes);
#endif //USE_XY
}

float Individual::EvaluateKinship(const Individual& other) {
	assert(chromosomeCount == other.chromosomeCount);


	int genomeSize=0;
	int totalKinship =0;
	for (uint i=0; i<chromosomeCount; i++) {

		totalKinship+=matChrom[i].EvaluateKinship(other.matChrom[i]);
		totalKinship+=matChrom[i].EvaluateKinship(other.patChrom[i]);
		totalKinship+=patChrom[i].EvaluateKinship(other.matChrom[i]);
		totalKinship+=patChrom[i].EvaluateKinship(other.patChrom[i]);
		genomeSize+=matChrom[i].LociCount();
	}

#ifdef USE_XY
	if (chromXY.chrom1) {
		totalKinship+=chromXY.chrom1->EvaluateKinship(other.chromXY.chrom1);
		totalKinship+=chromXY.chrom2->EvaluateKinship(other.chromXY.chrom2);
	}
#endif //USE_XY
	return (float)totalKinship/(float)(4*genomeSize);
}

void Individual::WritePedigreeMetaData(std::ostream& os) {
	static int statusConvertor[] = {1, 2, 0, 0};
	os<<pedID<<" "
		<<id<<" "
		<<matID<<" "
		<<patID<<" ";
	if (StandardPedigreeHeader) 
		os<<"0 0 0 "<<gender<<" 0 ";
	else
		os<<gender<<" ";
	os << statusConvertor[status] << " ";

	//Write each of the chromosomes to the line
	for (uint chID=0; chID<chromosomeCount; chID++) {

		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			//Let's make sure the minor allele frequencies are within the threshold
			if (MinorAlleleFreq(chID, locID) >= (minAlFreqThreshold - 0.0001)) {
				//If this is set as missing, return -1
				os << "0 0 ";
			}
		}
	}

#ifdef USE_XY
	if (chromXY.chrom1) {
		LocusManager<LocusXY>* loci = chromXY.chrom1->GetLoci();
		int locusCount = chromXY.chrom1->GetLocusCount();
		os<<" ";
		if (chromXY.chrom1) 
			//I need to add an XY version of this routine
			assert(0);
			chromXY.ShowGenotypesPedigree(os, minAlFreqThreshold, -1, ' ', NULL);
	}
#endif //USE_XY
	
  os << std::endl;
}

void Individual::WritePedigree(std::ostream& os, uint *genotypeCounts, bool writeOutcome) {
	static int statusConvertor[] = {1, 2, 0, 0};
	os<<pedID<<" ";								//Column #1
	os<<id<<" ";								//Column #2
	os<<matID<<" ";								//Column #3
	os<<patID<<" ";								//Column #4

	if (StandardPedigreeHeader) {
		os<<"0 ";			//uint nextOffspring;	//Column #5
		os<<"0 ";			//uint nextPatSib;		//Column #6
		os<<"0 ";			//uint nextMatSib;		//Column #7
	}
	os<<gender<<" ";			//uint sex;				//Column #8
	if (StandardPedigreeHeader) {
		os<<"0 ";			//uint probStat;		//Column #9
	}
	assert(status < 4);
	os << statusConvertor[status] << " ";
	if (writeOutcome)
		os<<outcome<<" ";
	int chromCount = chromosomeCount;
#ifdef USE_XY
	if (chromXY.chrom1) 
		chromCount--;
#endif
	int locPos = 0;
	//Write each of the chromosomes to the line
	for (int chID=0; chID<chromCount; chID++) {

		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			//Let's make sure the minor allele frequencies are within the threshold
			if (MinorAlleleFreq(chID, locID) >= (minAlFreqThreshold - 0.0001)) {
				//If this is set as missing, return -1
				if (missingData[chID][locID]) {
					if (genotypeCounts)
						genotypeCounts[locPos++ * 4]++;
					os << "0 0 ";
				}
				else {
					if (Individual::PhasedPedigrees) 
						os<<patChrom[chID][locID]+1<<" "<<matChrom[chID][locID]+1<<" ";
					else {
						int gt = GetGenotype(chID, locID);
						if (genotypeCounts)
							genotypeCounts[locPos++ * 4 + gt+1]++;
						switch(gt){
							case 0:
								os << "1 1 ";
								break;
							case 1:
								os << "1 2 ";
								break;
							case 2:
								os << "2 2 ";
								break;
						}
					}
				}
			}	
			else 
				cout<<"Skipping locus "<<locID + 1<<" frequency "<<MinorAlleleFreq(chID, locID)<<"\n";
		}
	}
#ifdef USE_XY
	if (chromXY.chrom1) {
		LocusManager<LocusXY>* loci = chromXY.chrom1->GetLoci();
		int locusCount = chromXY.chrom1->GetLocusCount();
		os<<" ";
		if (chromXY.chrom1) 
			chromXY.ShowGenotypesPedigree(os, minAlFreqThreshold, -1, ' ', genotypeCounts);
	}
#endif //USE_XY
	
  os << std::endl;
}



#ifdef USE_XY
int Individual::GetGenotypeXY(uint locus) {
	return chromXY.GetGenotype(locus);
}

bool Individual::ChangeGenotypeXY(uint locus, uint errorDir) {
	//This doens't make any sense with the current approach for the chromosomes
	cerr<<"Attempt to change genotype without verification to it's \"completeness\" has occured\n";
	exit(1);
}

bool Individual::ClearLocusXY(uint loc) {
	//This doens't make any sense with the current approach for the chromosomes
	cerr<<"Attempt to change genotype without verification to it's \"completeness\" has occured\n";
	exit(1);

}


void Individual::SetXX(AlleleSource<LocusXY>* x1, AlleleSource<LocusXY>* x2, PAR_Region<LocusXY>* par) {
	this->par = par;
	chromXY.Init(x1, x2, x1->GetLoci());
	assert(VerifyChromosomeIs(chromXY.chrom1, true));
	assert(VerifyChromosomeIs(chromXY.chrom2, true));
	chromXY.IsXX(true);
}

void Individual::SetXY(AlleleSource<LocusXY>* x, AlleleSource<LocusXY>* y, PAR_Region<LocusXY>* par) {
	this->par = par;
	assert(VerifyChromosomeIs(x, true));
	assert(VerifyChromosomeIs(y, false));
	
	chromXY.Init(x, y, x->GetLoci());
	assert(VerifyChromosomeIs(chromXY.chrom1, true));
	assert(VerifyChromosomeIs(chromXY.chrom2, false));
	chromXY.IsXX(false);
}


bool Individual::VerifyChromosomeIs(bool isX) {
	//Always should be X
	VerifyChromosomeIs(chromXY.chrom1, true);
	VerifyChromosomeIs(chromXY.chrom2, isX);
}

bool Individual::VerifyChromosomeIs(AlleleSource<LocusXY>* chromosome, bool isX) {
	stringstream ss1, ss2;
	LocusManager<LocusXY> *loci = chromosome->GetLoci();
	int locCount = loci->LocusCount();
	bool allClear = true;

	for (int i=0; i<locCount; i++) {
		LocusXY *loc = loci->At(i);
		ss1<<chromosome->At(i)<<" ";
		if (isX) {
			if (loc->type == LocusXY::Y_Only || loc->type == LocusXY::Y_Homolog) {
				if (chromosome->At(i)) {
					ss2<<"X ";
					allClear = false;
				}
				else 
					ss2<<"  ";
			}
			else 
				ss2<<"  ";
		
		} else {
			if (loc->type == LocusXY::X_Only || loc->type == LocusXY::X_Homolog) {
				if (chromosome->At(i)) {
					ss2<<"X ";
					allClear = false;
				}
				else 
					ss2<<"  ";
			}
			else 
				ss2<<"  ";
		}

	}
	if (!allClear)
		cerr<<ss1.str()<<"\n"<<ss2.str()<<"\n";
	return allClear;
}


#endif //USE_XY

void Individual::WriteContinuous(std::ostream&os, uint *genotypeCounts, bool doWriteStatus /*=true*/) {

	if (doWriteStatus)
		os << outcome << " ";

	int locPos = 0;
	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			//Let's make sure the minor allele frequencies are within the threshold
			if (MinorAlleleFreq(chID, locID) >= (minAlFreqThreshold - 0.0001)) {
	
				if (missingData[chID][locID]) {
					if (genotypeCounts)
						genotypeCounts[locPos++ * 4]++;
					os <<" -1";
				}			//If this is set as missing, return -1
				else {
					int gt=GetGenotype(chID, locID);
					if (genotypeCounts)
						genotypeCounts[locPos * 4 + gt + 1]++;
					os <<" "<<gt;
				}
			}
		}
	}

#ifdef USE_XY
	if (chromXY.chrom1) {
		LocusManager<LocusXY>* loci = chromXY.chrom1->GetLoci();
		int locusCount = chromXY.chrom1->GetLocusCount();
		os<<" ";
		if (chromXY.chrom1) 
			chromXY.ShowGenotypes(os, minAlFreqThreshold, -1, ' ', genotypeCounts);
	}
#endif //USE_XY
  os << std::endl;
}

void Individual::WriteMDR(std::ostream&os, uint *genotypeCounts, bool doWriteStatus /*=true*/) {

	if (doWriteStatus)
		os << status << " ";

	int locPos = 0;
	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			//Let's make sure the minor allele frequencies are within the threshold
			if (MinorAlleleFreq(chID, locID) >= (minAlFreqThreshold - 0.0001)) {
	
				if (missingData[chID][locID]) {
					if (genotypeCounts)
						genotypeCounts[locPos++ * 4]++;
					os <<" -1";
				}			//If this is set as missing, return -1
				else {
					int gt=GetGenotype(chID, locID);
					if (genotypeCounts)
						genotypeCounts[locPos * 4 + gt + 1]++;
					os <<" "<<gt;
				}
			}
		}
	}

#ifdef USE_XY
	if (chromXY.chrom1) {
		LocusManager<LocusXY>* loci = chromXY.chrom1->GetLoci();
		int locusCount = chromXY.chrom1->GetLocusCount();
		os<<" ";
		if (chromXY.chrom1) 
			chromXY.ShowGenotypes(os, minAlFreqThreshold, -1, ' ', genotypeCounts);
	}
#endif //USE_XY
  os << std::endl;
}

void Individual::ReadBinaryMdr(std::ifstream *genotypes, bool readMissingData, bool doWriteStatus, vector<uint> &chrSize) {
	//if (doWriteStatus)
	genotypes->read((char *)&status, 4);
		
	for (uint chID=0; chID<chromosomeCount; chID++) 
		patChrom[chID].ReadBinary(*genotypes, chrSize[chID], NULL, true);
#ifdef USE_XY
	if (chromXY.chrom1) {
		chromXY.chrom1->ReadBinary(*genotypes, NULL);
	}
#endif //USE_XY
	for (uint chID=0; chID<chromosomeCount; chID++) 
		matChrom[chID].ReadBinary(*genotypes, chrSize[chID], NULL, true);
#ifdef USE_XY
	if (chromXY.chrom2) {
		chromXY.chrom2->ReadBinary(*genotypes, NULL);
	}
#endif //USE_XY
	
	if (readMissingData) {
		for (uint chID=0; chID<chromosomeCount; chID++) {
			uint totalBlocks = (chrSize[chID] / boost::dynamic_bitset<>::bits_per_block) + 1;
			dynamic_bitset<>::block_type raw[totalBlocks];
			genotypes->read((char *)raw[0], (sizeof(dynamic_bitset<>::block_type)*totalBlocks));
			missingData[chID] = dynamic_bitset<>(&raw[0], &raw[totalBlocks]);
			missingData[chID].resize(chrSize[chID]);
		}
#ifdef USE_XY
		uint totalBlocks = (chromXY.chrom1->GetLocusCount() / boost::dynamic_bitset<>::bits_per_block) + 1;
		dynamic_bitset<>::block_type raw[totalBlocks];
		genotypes->read((char *)raw[0], (sizeof(dynamic_bitset<>::block_type)*totalBlocks));
		chromXY.missingData = dynamic_bitset<>(&raw[0], &raw[totalBlocks]);
		chromXY.missingData.resize(chromXY.chrom1->GetLocusCount());
#endif //USE_XY
	}
	
}
void Individual::WriteContinuous(std::ostream& meta, std::ostream&genotypes, bool writeMissingData, bool doWriteStatus /*=true*/) {

	//For now, we don't have a binary convention to leave off the status
	//if (doWriteStatus)
		genotypes.write((const char*)&outcome, 4);
		
	for (uint chID=0; chID<chromosomeCount; chID++) 
		patChrom[chID].WriteBinary(genotypes);
#ifdef USE_XY
	if (chromXY.chrom1) 
		chromXY.chrom1->WriteBinary(genotypes);
#endif //USE_XY

	for (uint chID=0; chID<chromosomeCount; chID++) 
		matChrom[chID].WriteBinary(genotypes);
#ifdef USE_XY
	if (chromXY.chrom2) 
		chromXY.chrom2->WriteBinary(genotypes);
#endif //USE_XY
	
	if (writeMissingData) {
		for (uint chID=0; chID<chromosomeCount; chID++) {
			BitSetType &misses = missingData[chID];
			uint totalBlocks = (misses.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
			dynamic_bitset<>::block_type raw[totalBlocks];

			//Convert the bitset to raw data
			to_block_range(misses, raw);

			genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	
		}
	}
#ifdef USE_XY
		BitSetType &misses = chromXY.missingData;
		uint totalBlocks = (misses.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
		dynamic_bitset<>::block_type raw[totalBlocks];

		//Convert the bitset to raw data
		to_block_range(misses, raw);

		genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	
#endif //USE_XY
}

void Individual::WriteMDR(std::ostream& meta, std::ostream&genotypes, bool writeMissingData, bool doWriteStatus /*=true*/) {

	//For now, we don't have a binary convention to leave off the status
	//if (doWriteStatus)
		genotypes.write((const char*)&status, 4);
		
	for (uint chID=0; chID<chromosomeCount; chID++) 
		patChrom[chID].WriteBinary(genotypes);
#ifdef USE_XY
	if (chromXY.chrom1) 
		chromXY.chrom1->WriteBinary(genotypes);
#endif //USE_XY

	for (uint chID=0; chID<chromosomeCount; chID++) 
		matChrom[chID].WriteBinary(genotypes);
#ifdef USE_XY
	if (chromXY.chrom2) 
		chromXY.chrom2->WriteBinary(genotypes);
#endif //USE_XY
	
	if (writeMissingData) {
		for (uint chID=0; chID<chromosomeCount; chID++) {
			BitSetType &misses = missingData[chID];
			uint totalBlocks = (misses.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
			dynamic_bitset<>::block_type raw[totalBlocks];

			//Convert the bitset to raw data
			to_block_range(misses, raw);

			genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	
		}
	}
#ifdef USE_XY
		BitSetType &misses = chromXY.missingData;
		uint totalBlocks = (misses.size() / boost::dynamic_bitset<>::bits_per_block) + 1;
		dynamic_bitset<>::block_type raw[totalBlocks];

		//Convert the bitset to raw data
		to_block_range(misses, raw);

		genotypes.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	
#endif //USE_XY
}

void Individual::WriteChromosome(std::ostream& os, uint chrID) {
	patChrom[chrID].ShowGenotypes(100, '\t');
#ifdef USE_XY
	/************ I'm not sure about this, these above don't do low MAF...that might be wrong */
	if (chromXY.chrom1) 
		chromXY.chrom1->ShowGenotypes(os, 100, '\t');
#endif //USE_XY
	matChrom[chrID].ShowGenotypes(100, '\t');
#ifdef USE_XY
	if (chromXY.chrom2) 
		chromXY.chrom2->ShowGenotypes(os, 100, '\t');
#endif //USE_XY
}

/**
 * I'm not completely certain this is correct. I need to visit the haploview manual for a refresher
 */
void Individual::WritePhased(std::ostream& os, uint indID, bool useFreqThreshold) {
	std::stringstream s1, s2;

	s1<<indID<<" ";				//Column #1
	s2<<indID<<" ";				//Column #1
	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		for (uint locID=0; locID<lociCount; locID++) {
			if (MinorAlleleFreq(chID, locID) >= (minAlFreqThreshold - 0.0001)) {
				s1<<patChrom[chID].At(locID) + 1<<" ";
				s2<<matChrom[chID].At(locID) + 1<<" ";
			}
		}
	}
#ifdef USE_XY
	if (chromXY.chrom1) 
		s1<<chromXY.chrom1;
	if (chromXY.chrom2)
		s2<<chromXY.chrom2;
#endif //USE_XY

 	os << s1.str() << std::endl;
 	os << s2.str() << std::endl;

}



#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(IndividualTest);
IndividualTest::IndividualTest() { }
IndividualTest::~IndividualTest() { 
	vector<Individual *>::iterator itr = individuals.begin();
	vector<Individual *>::iterator end = individuals.end();
int i=0;
	while (itr != end) {
		delete *itr;
		i++;
		itr++;
	}

/*	vector<Chromosome *>::iterator citr = chromosomes.begin();
	vector<Chromosome *>::iterator cend = chromosomes.end();
	while (citr!=cend) {
		delete *citr;
		citr++;
	}*/
}
void IndividualTest::setUp() {
	RecDistType recombIndexLookup;
	int locusCount = 10000;
	float mapPosition = 0.0;
	for (int i=0; i<locusCount; i++) {
		Locus l(0.2, 1, i, i+10);
		mapPosition += 0.2;
		l.MapPosition(mapPosition);
		loci.push_back(l);
		recombIndexLookup.Insert(mapPosition, i);
	}
	int lambda = 10;
	Chromosome *founder1 = new Chromosome(&loci, lambda, &recombIndexLookup);
	Chromosome *founder2 = new Chromosome(&loci, lambda, &recombIndexLookup);
	chromosomes.push_back(founder1);
	chromosomes.push_back(founder2);

	Chromosome *founder3 = new Chromosome(&loci, lambda, &recombIndexLookup);
	Chromosome *founder4 = new Chromosome(&loci, lambda, &recombIndexLookup);
	chromosomes.push_back(founder3);
	chromosomes.push_back(founder4);

	Chromosome *ex1	= new Chromosome(&loci, lambda, &recombIndexLookup);
	Chromosome *ex2 = new Chromosome(&loci, lambda, &recombIndexLookup);
	Chromosome *ex3 = new Chromosome(&loci, lambda, &recombIndexLookup);
	Chromosome *ex4 = new Chromosome(&loci, lambda, &recombIndexLookup);
	chromosomes.push_back(ex1);
	chromosomes.push_back(ex2);
	chromosomes.push_back(ex3);
	chromosomes.push_back(ex4);
	
	for (int i=0; i<locusCount; i++) {
		founder1->SetValue(i, 0);
		founder2->SetValue(i, 0);
		founder3->SetValue(i, 1);
		founder4->SetValue(i, 1);
		ex1->SetValue(i, 0);
		ex2->SetValue(i, 1);
		ex3->SetValue(i, 0);
		ex4->SetValue(i, 1);
	}

	Individual *fndInd1 = new Individual(1, 1, 1);
	fndInd1->SetChromosomalData(0, *founder1, *founder2);
	fndInd1->InitAlleleSource();
	individuals.push_back(fndInd1);							//[0] - founder
	
	Individual *fndInd2 = new Individual(2, 1, 1);
	fndInd2->SetChromosomalData(0, *founder3, *founder4);
	fndInd2->InitAlleleSource();
	individuals.push_back(fndInd2);							//[1] - Founder

	Individual *child1 = fndInd2->DLCross(fndInd1, 10);
	individuals.push_back(child1);							//[2] - Child 1
	Individual *child2 = fndInd2->DLCross(fndInd1, 11);
	individuals.push_back(child2);							//[3] - Child 2

	Individual *spouse1 = new Individual(10000, 1, 1);
	spouse1->SetChromosomalData(0, *ex1, *ex2);
	spouse1->InitAlleleSource();
	individuals.push_back(spouse1);							//[4] - Spouse
	
	Individual *spouse2 = new Individual(10001, 1, 1);
	spouse2->SetChromosomalData(0, *ex1, *ex2);
	spouse2->InitAlleleSource();
	individuals.push_back(spouse2);							//[5] - Spouse



	Individual *grandchild1 = child1->DLCross(spouse1, 100);
	individuals.push_back(grandchild1);						//[6] - Child
	Individual *grandchild2 = child1->DLCross(spouse1, 101);
	individuals.push_back(grandchild2);						//[7] - Child
	Individual *grandchild3 = child2->DLCross(spouse2, 102);
	individuals.push_back(grandchild3);						//[8] - Child
	Individual *grandchild4 = child2->DLCross(spouse2, 103);
	individuals.push_back(grandchild4);						//[9] - Child

	Individual *spouse3 = new Individual(20001, 1, 1);
	spouse3->SetChromosomalData(0, *ex1, *ex2);
	spouse3->InitAlleleSource();
	individuals.push_back(spouse3);							//[10] - Spouse
	
	Individual *spouse4 = new Individual(20002, 1, 1);
	spouse4->SetChromosomalData(0, *ex1, *ex2);
	spouse4->InitAlleleSource();
	individuals.push_back(spouse4);							//[11] - Spouse

	Individual *greatgc1 = grandchild1->DLCross(spouse3, 200);	
	individuals.push_back(greatgc1);						//[12] - Spouse
	Individual *greatgc2 = grandchild3->DLCross(spouse4, 201);
	individuals.push_back(greatgc2);						//[13] = Spouse
	
}

void IndividualTest::tearDown() { 
	
}

void IndividualTest::TestKinship() {
	//Test kinship with self
	float locusCount = 1000.0;
	float kinship = (float)individuals[0]->EvaluateKinship(*individuals[0]);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (Self)", 0.5, kinship, 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (none)", 0.0, (float)individuals[0]->EvaluateKinship(*individuals[1]), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (parent child)", 0.25, (float)individuals[1]->EvaluateKinship(*individuals[2]), 0.075);

	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (sibling)", 0.25, (float)individuals[2]->EvaluateKinship(*individuals[3]), 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (gran/child)", 0.125, (float)individuals[1]->EvaluateKinship(*individuals[6]), 0.125);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (cousin)", 0.125, (float)individuals[6]->EvaluateKinship(*individuals[8]), 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (uncle/niece)", 0.125, individuals[2]->EvaluateKinship(*individuals[8]), 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (uncle/niece)", 0.125, individuals[8]->EvaluateKinship(*individuals[2]), 0.05);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (second cousins)", 0.015625, individuals[12]->EvaluateKinship(*individuals[13]), 0.008);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (great grand)", 0.0625, individuals[12]->EvaluateKinship(*individuals[1]), 0.036);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Kinship (First Cousins once removed)", 0.0625, individuals[6]->EvaluateKinship(*individuals[13]), 0.03125);
}

#endif



}
