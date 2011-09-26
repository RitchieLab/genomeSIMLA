//
// C++ Interface: snpverificationmethod
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESESNPVERIFICATIONMETHOD_H
#define ESESNPVERIFICATIONMETHOD_H

#include "utility/utility.h"
#include "snpaligned.h"

namespace Genetics {
namespace ValidationFunctors {


/**
@brief Base class for the verification method functors used to determine whether or not to pass a given snp into a recipient

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpVerificationMethod{
public:
    SnpVerificationMethod() : snpsSeen(0) {}
    virtual ~SnpVerificationMethod(){}
	/**
	 * @brief Returns the number of snps observed by the method
	 */
	virtual uint GetSnpCount();

	/**
	 * @brief Resets the observed snp count and anything else associated with processing
	 */
	virtual void ResetSnpCount();
	
	/**
	 * @brief Evaluate a model
	 * @param model The model being evaluated
	 * @param position The position within the array where this model was found
	 */
	virtual bool EvaluateModel(SnpAligned *model, uint position) = 0;

	/**
	 * @brief Called immediately after evaluation is performed. 
	 * Override this with functionality that is required for after an evaluation has been made
	 */	
	virtual void PostEvaluation() {}

	/**
	 * @brief This is called right before evaluation is performed (1 time). 
	 * Override this for setting up things that need to be done shortly before the run.
	 */ 
	virtual void PreEvaluation() {}
protected:
	uint snpsSeen;
};

inline
uint SnpVerificationMethod::GetSnpCount() {
	return snpsSeen;
}

inline
void SnpVerificationMethod::ResetSnpCount() {	
	snpsSeen=0;
}

/**
@brief This functor passes any valid snp passed to it (valid being a pointer that isn't NULL)
*/
class AnySnp : public SnpVerificationMethod {
public:
	AnySnp() {}
	~AnySnp() {}
	bool EvaluateModel(SnpAligned *snp, uint position);


};
inline
bool AnySnp::EvaluateModel(SnpAligned *snp, uint position) {
	snpsSeen++;
	return true;
}

/**
@brief Passes only snps that are perfect
*/
class IsPerfect: public SnpVerificationMethod {
public:
	IsPerfect(const CaseControlStatus &status):  bitMask(status.affected) {}
	~IsPerfect() {}
	/**
	 * @brief Evaluate a model for whether it is perfect or not
	 * @param model the model being observed
	 * @param position The snp's index into the repository
	 */  
	bool EvaluateModel(SnpAligned *model, uint position);
protected:
	BitSetType bitMask; 		///<The mask used to determine affected/unaffected status
};

inline
bool IsPerfect::EvaluateModel(SnpAligned *snp, uint position) {
	snpsSeen++;
	if (snp)
		return snp->IsPerfect(bitMask) != 0;
	else
		return false;
}


class IsActive: public SnpVerificationMethod {
public:
	IsActive(const BitSetType &activeMask);
	~IsActive() {}

	bool EvaluateModel(SnpAligned *snp, uint position);	
protected:
	BitSetType activeMask;						///<Used to identify which bits are associated with affected individuals
};


inline
IsActive::IsActive(const BitSetType &activeMask) : activeMask(activeMask) {}

inline
bool IsActive::EvaluateModel(SnpAligned *snp, uint position) {
	snpsSeen++;
	assert(activeMask.size() > position);
	return activeMask[position];
}


}//Validators

}//ESE
#endif
