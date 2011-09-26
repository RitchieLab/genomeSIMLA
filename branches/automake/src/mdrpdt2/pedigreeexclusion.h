//
// C++ Interface: pedigreeexclusion
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_VALIDATORPEDIGREEEXCLUSION_H
#define MDRPDT_VALIDATORPEDIGREEEXCLUSION_H
#include <string>
#include <set>
#include <vector>
#include "pedigreerepositorytraverser.h"

namespace MdrPDT {

namespace Validator {


/**
@brief Exclude 0 or more pedigrees from the analysis

	@author Eric torstenson
*/
class PedigreeExclusion : public PedigreeRepositoryTraverser{
public:
    PedigreeExclusion(std::vector<std::string>& pedigrees);
    ~PedigreeExclusion();

	/**
	 * @brief Check for pedigree's presence in the exclusion list, and turn them off
	 */
	void Evaluate(Pedigree* ped);
	
	int ExclusionCount() { return exclusionCount; }
	std::string ExcludedPedigrees() { return excludedPedigrees; }
protected:
	std::set<std::string> exclusions;
	int exclusionCount;
	std::string excludedPedigrees;
};

}

}

#endif
