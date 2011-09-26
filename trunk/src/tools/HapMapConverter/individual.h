//
// C++ Interface: individual
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TOOLSINDIVIDUAL_H
#define TOOLSINDIVIDUAL_H
#include <vector>
#include "locusconverter.h"
namespace Tools {

using namespace std;

/**
	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class Individual{
public:
    Individual(vector<LocusConverter*>* loci);
    ~Individual();

	bool Parse(istream& file);
	string GetLabel();
	void SetLabel(const char *label);
	void SetGender(int gender);
	
	string GetGenotype(int idx);

protected:
	vector<LocusConverter*>* loci;
	vector<int> genotype1;
	vector<int> genotype2;
	string label;
	int gender;
};

inline
Individual::Individual(vector<LocusConverter*>* loci) : loci(loci) { }

inline
Individual::~Individual() { }

inline
void Individual::SetGender(int gender) {
	assert(gender>0);
	this->gender=gender;
}

inline
void Individual::SetLabel(const char *label) {
	this->label = string(label);
}

inline
string Individual::GetLabel() { 
	return label;
}

inline
bool Individual::Parse(istream& file) {
	int genotypeCount = loci->size();

	if (file.eof()) 
		return false;
	
	int gt=0;
	int sum=0;								///<For sanity checking
	for (int i=0; i<genotypeCount; i++) {
		file>>gt;
		sum+=gt;
		genotype1.push_back(gt);
	}
	if (sum > genotypeCount || sum < 1)
		return false;

	sum=0;	
	for (int i=0; i<genotypeCount; i++) {
		file>>gt;
		sum+=gt;
		genotype2.push_back(gt);
	}
	cout<<".";
	cout.flush();
	if (sum > genotypeCount || sum < 1)
		return false;
	return true;

}


string Individual::GetGenotype(int idx) {
	assert((size_t)idx < loci->size() && (size_t)idx < genotype1.size());
	return (*loci)[idx]->GetGenotype(genotype1[idx], genotype2[idx]); 
}


}

#endif
