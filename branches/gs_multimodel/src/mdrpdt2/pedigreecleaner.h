//
// C++ Interface: pedigreecleaner
//
// Description: 
// Detects genotype errors and cleans/removes various parts of the
// dataset according to user's specification
// 
// 
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_VALIDATORPEDIGREECLEANER_H
#define MDRPDT_VALIDATORPEDIGREECLEANER_H
#include <iostream>
#include <vector>
#include "familyrepository.h"
#include "utility/types.h"

namespace MdrPDT {

namespace Validator {

/**
@brief Scans pedigrees for various pedigree errrors and filters out genotypes, individuals and SNPS/genotypes based on user's suggestion. 

	@author Eric Torstenson
*/
class PedigreeCleaner{
public:
    PedigreeCleaner(int level, int threshold);
    ~PedigreeCleaner();

	/**
	 * @brief Evaluate over the contents of the repository and take action/generate report
	 */
	void Evaluate(PedigreeRepository *repo, std::ostream& os);

	/**
	 * @brief Evaluate a given individual, returning the number of bad loci encountered
	 */
	int Evaluate(Individual *ind, std::ostream& os, Utility::BitSetType& problemLoci);

	/**
	 * @brief Evaluate over a pedigree 
	 */
	int Evaluate(Pedigree* ped, std::ostream& os, std::ostream& summary);

protected:
	int level;
	int threshold;
	std::vector<int> errorCounts;

	void StripLocus(Pedigree *ped, int locus);
};

}

}

#endif
