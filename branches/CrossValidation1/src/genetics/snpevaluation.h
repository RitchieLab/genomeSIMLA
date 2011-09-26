//
// C++ Interface: snpevaluation
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONSNPEVALUATION_H
#define GENETICS_EVALUATIONSNPEVALUATION_H
#include "snpverificationmethod.h"
#include "permutationtest.h"

namespace Genetics {

namespace Evaluation {

using namespace Reporting;

/**
@brief Base class for SNP evaluation methods. 
This module simply evaluates each model in a repository based on the behavior of the method- but it doesn't attempt to populate a secondary repository. It is assumed that the method itself will have the various reports required. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpEvaluation : public Genetics::ValidationFunctors::SnpVerificationMethod {
public:
    SnpEvaluation(bool useModelThreshold = true) : distribution(NULL), useModelThreshold(useModelThreshold) {}
    ~SnpEvaluation() {}
	/**
	 * @brief Write the Evaluation results to the stream
	 */
	virtual void ReportResults( ostream *os)=0;

	/**
	 * @brief Set the Distribution object to be used
	 */
	virtual void SetDistribution(PermutationTestDist *dist);

	/**
	 * @brief Returns the distribution object used by the evaluation object
	 */
	virtual PermutationTestDist *GetDistribution();

	/**
	 * @brief Returns the top score over all runs
	 */
	virtual float GetTopScore(SnpAligned *&snp)=0;

  	/**
	 * @brief Set the general status for all individuals involved
	 */
	void SetOverallStatus(CaseControlStatus &stat);



protected:
	PermutationTestDist *distribution;					///<The distribution associated with the evaluation 
	CaseControlStatus overallStatus;					///<Used for calculating basic snp evaluation
	float overallRatio;									///<Used to set overall status ratio Aff/Unaffected
	bool useModelThreshold;								///<instructs the evaluation methods to calculate the threshold for each model
};


inline
void SnpEvaluation::SetDistribution(PermutationTestDist *dist) {
	distribution=dist;
}

inline
PermutationTestDist *SnpEvaluation::GetDistribution() {
	return distribution;
}

inline
void SnpEvaluation::SetOverallStatus(CaseControlStatus &stat) {
	overallStatus=stat;
	overallRatio = (float)stat.affected.count() / (float)stat.unaffected.count();
}

}

}

#endif
