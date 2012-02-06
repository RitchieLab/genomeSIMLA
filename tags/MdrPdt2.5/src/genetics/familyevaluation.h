//
// C++ Interface: familyevaluation
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFAMILYEVALUATION_H
#define GENETICS_EVALUATIONFAMILYEVALUATION_H
#include "utility/utility.h"
#include "familymember.h"

namespace Genetics {

namespace Evaluation {


/**
@brief Abstract base class describing functionality of an object that Evaluations a family
*/
class FamilyEvaluation {
public:
	FamilyEvaluation();
	virtual ~FamilyEvaluation();

	/**
	 * @brief Return the number of families were encountered during processing
	 */
	virtual uint GetFamilyCount();

	/**
	 * @brief Reset the count in case the same object is used more than once
	 */
	virtual void ResetFamilyCount();

	/**
	 * @brief Performs an evaluation to determine the usefulness of a given family member
	 */
	virtual bool Evaluate(FamilyMember *family, uint position) = 0;	
protected:
	uint familyCount;					///<Used to keep up with the number of families were visited
};

	
inline
FamilyEvaluation::FamilyEvaluation() : familyCount(0) {}

inline
FamilyEvaluation::~FamilyEvaluation() { }

inline
uint FamilyEvaluation::GetFamilyCount() {
	return familyCount;
}

inline
void FamilyEvaluation::ResetFamilyCount() {
	familyCount=0;
}



}

}

#endif
