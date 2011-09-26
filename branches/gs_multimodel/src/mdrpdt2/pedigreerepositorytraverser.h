//
// C++ Interface: pedigreerepositorytraverser
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_VALIDATORPEDIGREEREPOSITORYTRAVERSER_H
#define MDRPDT_VALIDATORPEDIGREEREPOSITORYTRAVERSER_H
#include "familyrepository.h"
#include "utility/types.h"
#include <utility>				///<std::pair

namespace MdrPDT {

namespace Validator {

/**
@brief Base class to help traverse pedigree repositories

	@author Eric Torstenson
*/
class PedigreeRepositoryTraverser {
public:
    PedigreeRepositoryTraverser();

    virtual ~PedigreeRepositoryTraverser();


	/**
	 * @brief Iterate through the contents of a given repository
	 */
	virtual void EvaluateRepository(PedigreeRepository *repo);

	/**
	 * @brief Perform tasks associated with a specific pedigree
	 * @note There is no need to call Evaluate(Individual), that is done by another function
	 */
	virtual void Evaluate(Pedigree* ped) { }

	/**
	 * @brief Perform tasks associated with a specific individual
	 */
	virtual void Evaluate(Pedigree* ped, Individual *ind) { }
	
	/**
	 * @brief Prepare to start an evaluation (basic book keeping)
	 */
	virtual void PrepEvaluation() { }

	/**
	 * @brief Finish up with the details of an evaluation
	 */
	virtual void PostEvaluation() { }
protected:
	/**
	 * @brief Interates over contents of a pedigree and calls the Evaluation for each individual
	 */
	virtual void _Evaluate(Pedigree* ped);

};





}

}

#endif
