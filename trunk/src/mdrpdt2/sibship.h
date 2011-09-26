//
// C++ Interface: sibship
//
// Description: 
//		Encapsulates the basic functionality of an extended sibship used in setting up
//		the genotype array for MDR-PDT analysis
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTSIBSHIP_H
#define MDRPDTSIBSHIP_H
#include <vector>
#include "utility/random.h"
#include "individual.h"

namespace MdrPDT {

/**
@brief Manages the production of valid DSPs based on sibships. 

When an individual is added, it will be added to the correct vector based on status. 
I

	@author 
*/
class Sibship{
public:
    Sibship();

	Sibship(const Sibship& other);

    ~Sibship();

	/**
	 * @add an individual to the sibship
	 */
	void AddIndividual(Individual *child);

	/**
	 * @brief return the number of DSPs associated with the sibship
	 * @note this does include those created by non-transmitted siblings, if there are any
	 */
	int GetDspCount();
	
	/**
	 * @brief populates the array with all of the DSP based genotypes and relevant meta data.
	 * @param start This is where to start writing data. 
	 * @param gen Random number generator, in case we are performing a PTest 
	 * @param xvSlice This is the cv slice in which the sibship belongs. This number is simply written into the header of the data
	 * @param permute This indicates wether or not the status should be permuted. 
	 * @param stride This is the distance between two neighboring snps for a given individual
	 * @return The pointer where the next sibship would start. This can be checked for correct
	 * stride by the calling routine.
	 */
	int WriteGenotypeData(char *start, int *meta, char *folds, int offset, Utility::Random& gen, int xvSlice, int stride, bool permute = false);
	char *WriteGenotypes(char *data, int *meta, char *folds, int offset, Individual *aff, Individual *unaff, int xvSlice, int stride);
	char *WriteGenotypes(char *data, int *meta, char *folds, int offset, Individual *affected, int xvSlice, int stride);
	char *WriteRandomizedTrio(char *data, int *meta, char *folds, int offset, Individual *affected, int xvSlice, int stride);
	/**
	 * @brief Returns the ids associated with the mother and father of the first child (which should be the same for all)
	 * @return t/f inicating that both parents have valid IDs
	 */
	bool GetParentIDs(std::string &father, std::string &mother);

	/**
	 * @brief Initialize the virtual siblings associated with each of the members of the sibship
	 */
	void InitVirtualSibs(std::ostream& os);

	/**
	 * @brief build the DSP count
	 */
	int CountDSPs();
protected:
	/**
	 * @brief Rebuild list of affecteds/unaffecteds based on random draw
	 */
	void Permute(std::vector<Individual*>& aff, std::vector<Individual*>& unaff, Utility::Random& gen);

	/**
	 * Members of these vectors are managed elsewhere....no need for managing locally
	 */
	std::vector<Individual*> affectedMembers;
	std::vector<Individual*> unaffectedMembers;

	/**
	 * @brief Cached DSP count based on real DSPs and those formed by virtual sibs.
	 * Is -1 until the count is performed.
	 */
	int effectiveDspCount;					///<Cached DSP count, based on real DSPs and those formed by virtual siblings	
};

}

#endif
