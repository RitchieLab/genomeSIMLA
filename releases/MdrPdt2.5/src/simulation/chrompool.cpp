//Chrompol.cpp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <iomanip>
#include <math.h>
#include "chrompool.h"
#include "utility/random.h"

namespace Simulation {

float ChromPool::defFre1 = 0.5;
float ChromPool::defFre2 = 0.5;
//uint ChromPool::errDirect = 0;
float ChromPool::errRate = 0.01;
bool ChromPool::randomizeAlleleFreq = true;

DistMappingFn *ChromPool::mappingFn = NULL;
//bool doUseMapFunctions	 = true;

using namespace std;
using namespace Utility;


//
// Default constructor
//
ChromPool::ChromPool(uint chromID, uint blockCount, float minR, float maxR) : 
	chromID(chromID), generationCount(0), lociCount(0), expressionCount(0), minRecomb(minR), maxRecomb(maxR), blockCount(blockCount), locSource("") {
	frequenciesDetermined = false;
	positionsDefined = false;
}

ChromPool::ChromPool(uint chromID, const char *locFilename) : chromID(chromID), generationCount(0), 
	lociCount(0), expressionCount(0), minRecomb(0), maxRecomb(0), blockCount(0), locSource(locFilename) {
	frequenciesDetermined = false;
	positionsDefined = false;
}

bool ChromPool::UseMapInfoFile() {
	return positionsDefined || mappingFn;
}

void ChromPool::BuildPool(uint expressionCount) {
	//If this fails, then we didn't initialize the loci vector
	cout<<"Building pool for chromosome: "<<chromID+1<<"\n";
	assert(loci.size() > 0);
	for (uint i=0; i<expressionCount; i++) {
		Chromosome ch(&loci);
		ch.InitLoci();
		pool.push_back(ch);
	}	

}	

/**
 * @brief Add a potential block of loci the chromosome
 * @param minCount The minimum number of snps 
 * @param maxCount The the maximum number of snps
 * @param minRecomb The minimum recombination Fraction between one locus and and the next one
 * @param maxRecomb The max recombination fraction between one locus and the next one
 * @param blockStrength Adjusts the block away from 0.5 chance. The larger the value, the weaker the "block" will be
 */
uint ChromPool::DefineBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax, float frequency) {
	assert(minSnpCount <= maxSnpCount);
	blockPrototypes.push_back(BlockDefinition(minSnpCount, maxSnpCount, blckMin, blckMax, snpMin, snpMax, frequency, blockPrototypes.size() + 1));
	return blockPrototypes.size();
}

void ChromPool::DefineDefaultBlock(ChromPool::BlockDefinition& block) {
	defaultBlock = block;
}

void ChromPool::DefineDefaultBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax) {
	assert(minSnpCount <= maxSnpCount);
	defaultBlock = ChromPool::BlockDefinition(minSnpCount, maxSnpCount, blckMin, blckMax, snpMin, snpMax, 100.0, 0);
}

ChromPool::BlockDefinition &ChromPool::DrawBlockDefinition() {
	float draw = Utility::Random::globalGenerator.drand();
	bool blockFound = false;
	uint blockCount = blockPrototypes.size();
	BlockDefinition &theBlock = defaultBlock;
	float currValue=0.0;
	for (uint i=0; i<blockCount && !blockFound; i++) {
		currValue += blockPrototypes[i].frequency;
		if (draw <= currValue) {
			theBlock = blockPrototypes[i];
			blockFound = true;
		}
	}
	return theBlock;		
}

double ChromPool::GetRecombinationFraction(uint locus) {
	double rec = 0.0;
	if (locus<lociCount)
		rec = loci[locus].RecombinationFraction();
	return rec;
}


void ChromPool::WriteMarkerInfo(std::ostream& os, DistMappingFn *f) {
	WriteMarkerInfo(os, 0, lociCount, f);
}

void ChromPool::WriteMarkerInfo(std::ostream& os, uint first, 
		uint last, DistMappingFn *f) {

	if (last == 0 || last > lociCount)
		last = lociCount;

	for(uint i=first; i<last; i++) 
		loci[i].WriteMarkerInfo(os);

}

uint ChromPool::InitLoci() {
	if (locSource.length() == 0) {
		
		int lPos = 0;					///<The physical location of the locus
		int markerID = 0;
		//ssert(minRecomb<=maxRecomb);
		cout<<"Initializing loci for chromosome: "<<chromID<<"\n";
		for (uint i=0; i<blockCount; i++){ 
			BlockDefinition &block = DrawBlockDefinition();
			uint snpCount = block.GetSnpCount();
			float recomb = block.GetBlockRecombination();
			//float rDiff = maxRecomb - minRecomb;
	
			stringstream blckID;
			blckID<<"Blck"<<i<<"x"<<block.id;
	
			for (uint snpID=0; snpID<snpCount; snpID++) {
				lPos = 0;
				if (mappingFn) 
					lPos = mappingFn->GetMapPosition(recomb);

				Locus l(recomb, errRate, ++markerID, lPos);
				if (snpID == 0) 
					l.SetDescription(blckID.str().c_str());
				
				//Determine how we assign allele frequencies
				if (randomizeAlleleFreq)
					l.RandomizeFreq(defFre1, defFre2);
				else
					l.AssignFreq(defFre1, defFre2);
		
				//Add the loci to the collection
				loci.push_back(l);
				//Adjust the recombination value for the next locus
				//recomb=minRecomb+(rDiff * Utility::Random::globalGenerator.drand());
				recomb = block.GetSnpRecombination();
				
			}
		}
		lociCount = loci.size();
	}
	else {
		cout<<"Attempting to load loci from the file: "<<locSource<<"\n";
		ifstream lf(locSource.c_str(), ios_base::in);
		if (lf.is_open()) {
			LoadLoci(lf);
		} 
		else {
			cout<<"Unable to open locus file\n";
			lociCount = 0;
		}
	}
	
	return lociCount;
}

uint ChromPool::CountColumns(const char *line) {
	stringstream ss(line);
	string vals;

	uint count = 0;

	while (!ss.eof()) {
		ss>>vals;	
		count++;
	}
	if (count > 0)
		count--;
	
	return count;
}

uint ChromPool::LoadLoci(istream& data) {
	frequenciesDetermined = true;

	lociCount = 0;
	
	char line[4096];
	data.getline(line, 4096);			///<The chromosome id
	data.getline(line, 4096);			///<The count of loci
	data.getline(line, 4096);			///<The header information

	string id;

	while (!data.eof()) {
		Locus loc;
		//data>>id;
		if (!data.eof()) {
			data>>loc;

			if (loc.Valid())
				loci.push_back(loc);

			//For debugging
			//cout<<id<<" "<<loc<<" \n";
			
		}
	}
	lociCount = loci.size();

	//Well, let's check the last locus for it's position to determine if there are valid positions 
	if (lociCount > 0)
		positionsDefined = loci[lociCount - 1].GetLocation() > 0;

	return lociCount;



}



void ChromPool::LoadPhased(istream& data, uint currGen) {
	uint expectedLineLength = ((lociCount * 2) + 2) * 2;
	char line[expectedLineLength];
	data.getline(line, expectedLineLength);
	data.seekg(0, ios::beg);

	uint columnCount = CountColumns(line);
	if (columnCount != lociCount + 2) {
		cout<<"Grrr! "<<columnCount<<" "<<lociCount + 2<<"\n";
		throw ("wrong number of SNPs!");
		//cout<<"The current gene pool doesn't match the configuration. Please check correct the settings or rerun the initial population once again.\n\tLoci Count in File: "<<((float)(columnCount-6)/2.0)<<" - Configuration file specifies: "<<lociCount<<"\n";
	}

	while (!data.eof()) {
		string var;
		
		//Ignore the first 6 columns
		for (uint i=0; i<2; i++) 
			data>>var;

		if (!data.eof() ) {

			//Get the two strands of DNA 
			Chromosome c(&loci);

			for (uint i=0; i<lociCount; i++) {
				int strand1;
				data>>strand1;
				c.SetValue( i, strand1 );
			}
			pool.push_back(c);
		}
	}
	
	lociCount = loci.size();

	cout<<"Previous gene pool has been loaded. Configuration parameters might not explicitly describe the current set of individuals. \n";
}

/**
 * Performs generational mating/replacement
 * @param genCount The number of generations to move forward. 
 * @return The current generation after we are finished
 */
uint ChromPool::AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f) {
	int firstChromIndex, secondChromIndex;
	frequenciesDetermined = false;
	
	cout<<"\nAdvancing chromosome "<<chromID<<" ("<<(pool.size()/2)<<") : ";
	for(uint currGen=0; currGen < genCount; currGen++){   
		PoolType lastpool = pool;
	 	int numChroms = lastpool.size();

		pool.clear();
		uint endSize = (*f)(++generationCount);

		for(uint currChrom=0; currChrom < endSize; currChrom++){
			firstChromIndex = int(Random::globalGenerator.drand() * numChroms);
			secondChromIndex = int(Random::globalGenerator.drand() * numChroms);
			pool.push_back(lastpool[firstChromIndex].Cross(lastpool[secondChromIndex]));
    	}
		cout<<"*";
		cout.flush();
  	}
	cout<<" ("<<(pool.size()/2)<<")\n";
	usedIndividuals.clear();
	return generationCount;
}




bool ChromPool::ForceAlleleFrequency(uint locusID, float al1, float al2) {
	bool success = locusID < lociCount;

	if (success) {
		cout<<"Setting allele frequency for chromosome: "<<chromID<<"\n";
		loci[locusID].AssignFreq(al1, al2);
	}
	else {
		cout<<"Unable to set frequency for chromosome "<<chromID<<"\n";
		cout<<"There is no locus "<<locusID<<" (Total loci "<<lociCount<<")\n";
	}

	return success;
}



/**
 * 0 based
 */
bool ChromPool::GetAlleleFrequencies(uint locus, float &af1, float &af2) {
	if (locus<lociCount) {
		af1 = loci[locus].Freq1();
		af2 = loci[locus].Freq2();
	}
	
	return locus<lociCount;
}

void ChromPool::CalculateAlleleFrequencies() {
	if (frequenciesDetermined)
		return;

	uint indCount = pool.size();
	//locix2 array for each allele
	uint freq[lociCount][2];
	for (uint i=0; i<lociCount; i++) {
		freq[i][0]=0;
		freq[i][1]=0;
	}

	//os << setw(width) << right << "Loc  "<< setw(width) << " All1 "<< setw(width) << "All2 "<< setw(width) << "Recomb   " << std::endl;
	for (uint ind = 0; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];
		for(uint i=0; i<lociCount; i++)
			freq[i][curChrom.At(i)]++;
	}

	for (uint i=0;i<lociCount; i++) {
		double al1 = (double)(freq[i][0])/(double)(indCount);
		double al2 = (double)(freq[i][1])/(double)(indCount);;
		loci[i].AssignFreq(al1, al2);
	}

	frequenciesDetermined = true;
	
}

void ChromPool::SaveLoci(ostream& os) {
    ios_base::fmtflags old_settings = os.flags();
	os.setf(ios::fixed, ios::floatfield);
	
	int width=16;

	os << "Locus Log for Chromosome "<<chromID<<"\n";
	os << lociCount<<" Loci\n";

	if (lociCount > 0) {
		//os << setw(width)<<"Locus ID" << " ";
		loci[0].WriteHeader(os, width);
	}
	for(uint i=0; i<lociCount; i++){
		//os << setw(width) << i+1 << " ";
		loci[i].WriteFormatted(os, width);
	}
  	os.setf(old_settings); 
}



//Let's dump the gene pool to a stream. We'll use haploview format
uint ChromPool::SaveAsPhased(std::ostream& os, uint first, uint count) {
	if (first == pool.size())
		return 0;

	PoolType::iterator i;
	PoolType::iterator end=pool.end();

	uint last=0;
	if (count == 0)
		last = pool.size() - first;
	else {
		last = first+count;
		if (last > pool.size())
			last = pool.size();
	}


	uint id=first;
	for (i=pool.begin(); i!=end && i+1 != end; ) {
		os<<++id<<" 1 ";
		(*i++).WritePedFormat(os, first, last);

		os<<id<<" 1 ";
		(*i++).WritePedFormat(os, first, last);
	}

	return last;
}

void ChromPool::SaveAsPhased(std::ostream& os) {
	PoolType::iterator i;
	PoolType::iterator end=pool.end();

	uint id=0;

	for (i=pool.begin(); i!=end && i+1 != end; ) {
		//i is pedigree. 1 is individual ID. All the rest of metadata is 0
		os<<++id<<" 1 ";
		(*i++).WritePedFormat(os, 0, 0);

		os<<id<<" 1 ";
		(*i++).WritePedFormat(os, 0, 0);
	}
}

float ChromPool::GetErrorRate(uint locus) {
	assert(locus<lociCount);

	return loci[locus].ErrorRate();
}



void ChromPool::ReturnIndividual( Individual *ind) {
	IndLookupType::iterator pos = usedIndividuals.find(ind->GetPoolID(chromID));
	IndLookupType::iterator end = usedIndividuals.end();

	if (pos != end)
		usedIndividuals.erase(pos);
	
}

void ChromPool::DrawIndividual(Individual* ind) {
	uint ChromCount = pool.size();
	uint chromID1;
	uint chromID2;

	char id[128];

	//Make sure we don't redraw any individual
	IndLookupType::iterator pos;
	IndLookupType::iterator end = usedIndividuals.end();

	uint attempts = 0;

	do {
		chromID1 = int(Random::globalGenerator.drand() * ChromCount);
		chromID2 = int(Random::globalGenerator.drand() * ChromCount);
		sprintf(id, "%dx%d", chromID1, chromID2);
		pos = usedIndividuals.find(id);
		if (attempts++>5000) {
			cout<<"We are having trouble getting individual "<<ind->GetID()<<" from chromosome: "<<chromID<<". \nThe list of used individuals is "<<usedIndividuals.size()<<" long\n";
			cout<<"Go back to your data and make sure that your pool is big enough to support the number of datasets and individuals that are being requested. With your current settings, the pool is: "<<ChromCount<<"\n";
			abort();
		}
	} while (pos != end);

	ind->SetPoolID(chromID, id);
	ind->SetChromosomalData(chromID, pool[chromID1], pool[chromID2]);

	//usedIndividuals[id]=ind->GetID();			//or something like this
	usedIndividuals[id]=1;
}



// Clears the chromosomal pool and individuals in this population
// Arg: none
// Ret: none
void ChromPool::Clear(){
	usedIndividuals.clear();
	pool.clear();
}




}
