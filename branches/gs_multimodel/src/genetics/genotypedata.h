//
// C++ Interface: snp
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSSNPINDALIGNED_H
#define GENETICSSNPINDALIGNED_H
#include <iostream>
#include "utility/utility.h"

namespace Genetics {
using namespace Utility;
using namespace std;
/**
@brief This is the functional element for storing genetic information for a family member. 
We first load data into these structures and build SnpAligned structures once we figure out 
what exactly we need.


	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenotypeData{
public:
    GenotypeData(TypeLookup *parser);

    ~GenotypeData();

	/**
	 * @brief Count the number of individuals associated with a given genotype
	 */
	uint CountGenotypes(int genotype);

	/**
	 * @brief Return the total number of individuals for this and other genotypes
 	 * @todo This is a terrible name for this function! Rename it
	 */
	uint GetGenotypeNumber(int genotype);
	/**
	 * @brief Returns the label for the genotype at position, position
	 */  
	string GetGenotypeLabel(uint position);
	
	/**
	 * @brief Returns the index of the genotype for position, position
	 */
	int GetGenotypeIndex(uint position);
	
	/**
	 * @brief Turns on the individual for a given string
	 */
	void SetGenotype(uint position, const char *encodedGenotype);
	void SetGenotype(uint position, uint idx);
	/**
	 * @brief Ensure that the storage has adequate storage for the current value
	 */
	void VerifyDimensions(int idx, uint position);

	/**
	 * @brief Generate a descriptive report of the genetic data for the individual
	 */
	ostream *Report(ostream *os);

	/**
	 * @brief This strips away data from the genotypes where data is missing in "other". 
	 * This is used in the generation of DSPs
	*/
	GenotypeData *BuildMaskedGT(GenotypeData *other);

	/**
	 * @brief Returns the number of genotypes are present in the structure
	 */
	uint CountGenotypesPresent();

	void ZeroGenotype(uint position);

protected:
	/**
     * @brief Contains a bitset for each of the genotypes (possible variables)
	 * For each genotype/discrete variable, there will be a bitset whose values represent individuals who
	 * have this genotype at a given SNP. 
	 */
	GenotypeArray genotypePresent;

	/**
	 * @brief This is used to lookup the various indices for a given genotype 
	 */
	TypeLookup *genoLookup;	
};

inline
uint GenotypeData::CountGenotypesPresent() {
	if (genotypePresent.size() == 0)
		return 0;
	return genotypePresent[0].individuals.count();
	//return genotypePresent[0].individuals.size() - genotypePresent[0].individuals.count();
} 

inline
GenotypeData::GenotypeData(TypeLookup *lkup) : genoLookup(lkup) { }

inline
GenotypeData::~GenotypeData() {} 



inline
void GenotypeData::VerifyDimensions(int idx, uint position) {

	if ((uint)idx > genotypePresent.size())						//Switched to < from <= in order to avoid resizing over and over
		genotypePresent.resize(idx);

	uint count=genotypePresent.size();
	
	if (genotypePresent[0].individuals.size() < position)
		genotypePresent[0].individuals.resize(position, true);

	for (uint i=1; i<count; i++) 
		if (genotypePresent[i].individuals.size() < position)
			genotypePresent[i].individuals.resize(position);
}

inline
ostream *GenotypeData::Report(ostream *os) {
	if (genotypePresent.size() == 0)	{
		*os<<"no genotype data present\n";
		return os;
	}
	uint count=genotypePresent[0].individuals.size();

	for (uint i=0; i<count; i++) 
		*os<<GetGenotypeLabel(i)<<"  ";

	return os;
}

inline
uint GenotypeData::GetGenotypeNumber(int genotype) {
	if (genotypePresent.size() == 0)
		return 0;
	return genotypePresent[genotype].individuals.size();
}

inline
uint GenotypeData::CountGenotypes(int genotype) {
	assert((uint)genotype < genotypePresent.size());
	return genotypePresent[genotype].individuals.count();
}

inline
string GenotypeData::GetGenotypeLabel(uint position) {
	int idx = -1;

	uint value=0;
	for (; idx == -1 && value < genotypePresent.size(); value++) 
		if (genotypePresent[value].individuals[position])
			idx=value;
	if (idx==-1)
		idx=0;
	return genoLookup->GetValue(idx);
}

	
	
inline
void GenotypeData::SetGenotype(uint position, uint idx) {
	VerifyDimensions(idx+1, position+1);
	genotypePresent[idx].individuals[position]=idx != 0;
}

inline
void GenotypeData::SetGenotype(uint position, const char *encodedGenotype) {
	int idx=genoLookup->GetValue( encodedGenotype);
	VerifyDimensions(idx+1, position+1);
	genotypePresent[idx].individuals[position]=idx != 0;
}


}

#endif
