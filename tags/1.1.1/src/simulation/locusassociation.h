//
// C++ Interface: locusassociation
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLOCUSASSOCIATION_H
#define SIMULATIONLOCUSASSOCIATION_H
#include "utility/types.h"
#include "utility/random.h"
#include <vector>

namespace Simulation {

/**
	@brief Aides in the production of associations at the initialization of the chromosomes
	@Note The class, LocusAssoc represents a single (or multilocus) association
	@author 
*/

class LocusAssociation{
public:
	/**
	 * @brief Construction
	 * @param Position The locus position (index) Should be greater than 0
	 * @param Level is the number of previous SNPs depended on by the local SNP
	 */
    LocusAssociation(int position, int level);
    ~LocusAssociation();

	/**
	 * Returns 0/1 based on the precedent and a random draw from r
	 */
	bool GetAllele(Utility::Random& r);

	bool Load(std::istream& input);		///<Load the appropriate numbers of values from the input file
	Utility::BitSetType precedent;			///<Clients use this to set the preceding alleles.
	
	float GetFrequencyZero();		///<Returns the percentage of 0s observed;
	void ResetObservations();		///<In case we need to reset the observations for frequencies
protected:	
	int position;					///<locus Index
	int level;						///<Number of loci used to determine the allele frequency
	std::vector<float> probability;		///<Probability of a 0-The size of probability is level^2
	int observedZeros;				///<Count the number of 0s observed
	int totalObservations;			///<Total number of observations
};

class LocusAssociationGrid {
public:
	void ResetObservations();
	int Load(const char *filename);
	Utility::BitSetType GenerateChromosome(Utility::Random& rnd);
protected:
	std::vector<LocusAssociation> loci;
	int gridDepth;
};

}
#endif
