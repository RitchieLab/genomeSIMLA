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

namespace Simulation {

using namespace std;
using namespace Utility;

bool Individual::StandardPedigreeHeader = true;

Individual::Individual(uint id, uint pedID, uint chrCount /*=1*/) : id(id), pedID(pedID), patID(0), matID(0), includeInDataset(true), status(0), chromosomeCount(chrCount) {
	matChrom = new Chromosome[chrCount];
	patChrom = new Chromosome[chrCount];
	genePoolIDs = new string[chrCount];
	missingData = new BitSetType[chrCount];
}

Individual::Individual() : id(0), pedID(0), patID(0), matID(0), includeInDataset(true), status(0), matChrom(NULL), patChrom(NULL), chromosomeCount(0), genePoolIDs(NULL) { }

void Individual::Init(uint id, uint pedID, uint chrCount /*=1*/)  {
	id = id;
	pedID = pedID;
	chromosomeCount = chrCount;
	matChrom = new Chromosome[chrCount];
	patChrom = new Chromosome[chrCount];
	if (genePoolIDs)
		delete[] genePoolIDs;
	genePoolIDs = new string[chrCount];
}

bool Individual::DoIncludeInDataset() {
	return includeInDataset;
}

void Individual::DoIncludeInDataset( bool doinclude){
	includeInDataset=doinclude;
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

void Individual::SetPedigreeMeta(uint patID, uint matID) {
	this->patID=patID;
	this->matID=matID;
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
	matChrom = m;
	patChrom = p;
	for (uint i=0; i<chromosomeCount; i++)
		missingData[i].resize(p[i].LociCount(), false);
}


void Individual::SetChromosomalData(uint chromID, Chromosome &p, Chromosome &m) {
	assert(chromID < chromosomeCount);
	matChrom[chromID]=m;
	patChrom[chromID]=p;
	missingData[chromID].resize(p.LociCount(), false);
}


void Individual::SetStatus(bool isAffected) {
	status = isAffected;
}

bool Individual::IsAffected() {
	return status != 0;
}


bool Individual::ClearLocus( uint chrID, uint locus) {
	if (missingData[chrID][locus]) 
		return false;
	else 
		missingData[chrID][locus] = true;
	return true;
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

bool Individual::ApplyStatus(PenetranceModel *model) {
	uint modelSize = model->GetModelSize();
	vector<uint> genotypes;
	for (uint i=0; i<modelSize; i++)  {
		genotypes.push_back(GetGenotype(model->GetChromID(i), model->GetLocus(i)));
	}

	bool isAffected = model->IsAffected(genotypes);
	
	SetStatus(isAffected);

	return isAffected;	
}

Individual *Individual::Cross(Individual *father, uint id) {
	Chromosome *pDNA = father->Cross();
	Chromosome *mDNA = Cross();
	Individual *child = new Individual(id, father->GetPedigreeID(), chromosomeCount);
	child->SetPedigreeMeta(father->GetID(), GetID());
	child->SetChromosomalData(pDNA, mDNA);
	return child;
}

Chromosome *Individual::Cross() {
	Chromosome *newChr = new Chromosome[chromosomeCount];
	for (uint i=0; i<chromosomeCount; i++)
		newChr[i] = patChrom[i].Cross(matChrom[i]);

	return newChr;
}

int Individual::GetGenotype(uint chrID, uint locus) {

	return matChrom[chrID][locus] + patChrom[chrID][locus];
}

int Individual::GetSnpCount() {
	uint lociCount = 0;
	for (uint chID=0; chID<chromosomeCount; chID++) 
		lociCount+=patChrom[chID].LociCount();

	return lociCount;
}

void Individual::WritePedigree(std::ostream& os, uint *genotypeCounts) {
	static int statusConvertor[] = {1, 2, 0};
	os<<pedID<<" ";								//Column #1
	os<<id<<" ";								//Column #2
	os<<matID<<" ";								//Column #3
	os<<patID<<" ";								//Column #4

	if (StandardPedigreeHeader) {
		os<<"0 ";			//uint nextOffspring;	//Column #5
		os<<"0 ";			//uint nextPatSib;		//Column #6
		os<<"0 ";			//uint nextMatSib;		//Column #7
	}
	if (id<3)
		os<<id<<" ";
	else
		os<<"0 ";			//uint sex;				//Column #8
	if (StandardPedigreeHeader) {
		os<<"0 ";			//uint probStat;		//Column #9
	}
	os << statusConvertor[status] << " ";

	int locPos = 0;
	//Write each of the chromosomes to the line
	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			//If this is set as missing, return -1
			if (missingData[chID][locID]) {
				if (genotypeCounts)
					genotypeCounts[locPos++ * 4]++;
				os << "0 0 ";
			}
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
	
  os << std::endl;
}

void Individual::WriteMDR(std::ostream&os, uint *genotypeCounts) {
	os << status << " ";

	int locPos = 0;
	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {

			if (missingData[chID][locID]) {
				if (genotypeCounts)
					genotypeCounts[locPos++ * 4]++;
				os << "-1";
			}			//If this is set as missing, return -1
			else {
				int gt=GetGenotype(chID, locID);
				if (genotypeCounts)
					genotypeCounts[locPos * 4 + gt + 1]++;
				os <<" "<<gt;
			}
		}
	}
  os << std::endl;
}

/**
 * I'm not completely certain this is correct. I need to visit the haploview manual for a refresher
 */
void Individual::WritePhased(std::ostream& os, uint indID) {
	std::stringstream s1, s2;

	s1<<indID<<" ";				//Column #1
	s2<<indID<<" ";				//Column #1

	for (uint chID=0; chID<chromosomeCount; chID++) {
		uint lociCount=patChrom[chID].LociCount();
		os << " ";
		for (uint locID=0; locID<lociCount; locID++) {
			s1<<patChrom[chID].At(locID) + 1<<" ";
			s2<<matChrom[chID].At(locID) + 1<<" ";
		}
	}
	
 	os << s1.str() << std::endl;
 	os << s2.str() << std::endl;

}


}
