//
// C++ Interface: minimalchromosome
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONMINIMALCHROMOSOME_H
#define SIMULATIONMINIMALCHROMOSOME_H

#include <vector>
#include <map>
#include "utility/types.h"
#include "utility/random.h"
#include "chromosome.h"

namespace Simulation {

using namespace std;

/**
To be used when the pool has been "closed". 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

class MinChrom {
public:
	struct Common {
		Utility::Random *generator;
		LocusArray *underlyingLoci;					///<Describes various loci
		double poissonLambda;						///<Lambda to be used by poisson distribution
		RecDistType *recombIndexLookup;				///<Quickly lookup the index at a given location
	};	
	
	MinChrom(vector<uint> *loci, Common *details);
	MinChrom(const MinChrom& other);
	
	MinChrom &operator=(const MinChrom& other);
	MinChrom &operator=(Chromosome &other);

	int At(uint locIdx);
	int operator[] (uint locusIndex);

	Common *GetCommonData() { return commonData; }

	uint GetModelIndex(uint i);
	uint GetModelSize();
protected:
	map<uint, bool> chromosomalData;				///<Locus-Index->allele 0/1
	vector<uint> *modelLoci;						///<The indices associated with the model from this chromosome
	Common *commonData;								///<To avoid keeping the same info multiple times

};

inline
MinChrom::MinChrom(vector<uint> *loci, Common *details) : modelLoci(loci), commonData(details) { }

inline
MinChrom::MinChrom(const MinChrom &other) {
	modelLoci 		= other.modelLoci;
	commonData 		= other.commonData;
	chromosomalData = other.chromosomalData;
}

inline
MinChrom &MinChrom::operator=(const MinChrom& other) {
	modelLoci 		= other.modelLoci;
	commonData 		= other.commonData;
	chromosomalData = other.chromosomalData;
	return *this;
}

inline
int MinChrom::At(uint locusIdx) {
	map<uint, bool>::iterator itr = chromosomalData.find(locusIdx);

	if (itr == chromosomalData.end()) {
		cout<<"No data at: "<<locusIdx<<" when querying a minimal chromosome for data\n";
		assert(0);
		return -1;
	}
	return itr->second;
}

inline
int MinChrom::operator[](uint locusIdx) {
	return At(locusIdx);
}

inline
MinChrom &MinChrom::operator=(Chromosome &other) {
	//We'll assume that the common stuff has already been set
	uint modelSize = modelLoci->size();
	chromosomalData.clear();
	
	for (uint i=0; i<modelSize; i++) {
		int b = modelLoci->at(i);
		chromosomalData[b] = other.At(b);
	}
	return *this;
}

}

#endif
