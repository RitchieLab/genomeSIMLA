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
#include "ptestdistribution.h"
#include "genetics/snpaligned.h"
#include "modelstatistics.h"

namespace Genetics {

namespace Evaluation {

using namespace Distribution;

/**
@brief Base class for SNP evaluation methods. 
This module simply evaluates each model in a repository based on the behavior of the method- but it doesn't attempt to populate a secondary repository. It is assumed that the method itself will have the various reports required. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpEvaluation : public Genetics::ValidationFunctors::SnpVerificationMethod {
public:
	
    SnpEvaluation(bool useModelThreshold = true) : 
		fitnessDist(NULL), 
		orDist(NULL), 
		peDist(NULL),
		lrDist(NULL),
		useModelThreshold(useModelThreshold) {}

    ~SnpEvaluation() {}
	/**
	 * @brief Write the Evaluation results to the stream
	 */
	virtual void ReportResults(ostream *os)=0;

	/**
	 * @brief Set the Distribution object to be used
	 */
	virtual void SetFitnessDist(PTestDistribution *dist);
	virtual void SetOddsRatioDist(PTestDistribution *dist);
	virtual void SetPEDist(PTestDistribution *dist);
	virtual void SetLRDist(PTestDistribution *dist);

	/**
	 * @brief Returns the distribution object used by the evaluation object
	 */
	virtual PTestDistribution *GetFitnessDist();
	virtual PTestDistribution *GetOddsRatioDist();
	virtual PTestDistribution *GetPEDist();
	virtual PTestDistribution *GetLRDist();

	/**
	 * @brief Returns the top score over all runs
	 * @param snp The model will be put here
	 */
	virtual ModelStatistics GetModelStats(SnpAligned *&snp)=0;

	virtual ModelStatistics GetModelStats(uint modelSize, SnpAligned *&snp)=0;

	virtual FoldStatistics GetFoldStatistics(uint foldID, SnpAligned *snp)=0;

	/**
	 * @brief Creates the value used for distribution stuff
	 */
	virtual float GetAverageTesting(SnpAligned *&topSnp)=0;

  	/**
	 * @brief Set the general status for all individuals involved
	 */
	void SetOverallStatus(CaseControlStatus &stat);



protected:
	PTestDistribution *fitnessDist;						///<The distribution associated with the evaluation 
	PTestDistribution *orDist;							///<Distribution associated wtih the odds ratio
	PTestDistribution *peDist;
	PTestDistribution *lrDist;
	CaseControlStatus overallStatus;					///<Used for calculating basic snp evaluation
	float overallRatio;									///<Used to set overall status ratio Aff/Unaffected
	bool useModelThreshold;								///<instructs the evaluation methods to calculate the threshold for each model
};


inline
void SnpEvaluation::SetFitnessDist(PTestDistribution *dist) {
	fitnessDist=dist;
}

inline
void SnpEvaluation::SetOddsRatioDist(PTestDistribution *dist) {
	orDist = dist;
}

inline
void SnpEvaluation::SetPEDist(PTestDistribution *dist) {
	peDist = dist;
}

inline
void SnpEvaluation::SetLRDist( PTestDistribution *dist ){
	lrDist = dist;
}

inline
PTestDistribution *SnpEvaluation::GetFitnessDist() {
	return fitnessDist;
}

inline
PTestDistribution *SnpEvaluation::GetLRDist() {
	return lrDist;
}

inline
PTestDistribution *SnpEvaluation::GetOddsRatioDist() {
	return orDist;
}

inline
PTestDistribution *SnpEvaluation::GetPEDist() {
	return peDist;
}

inline
void SnpEvaluation::SetOverallStatus(CaseControlStatus &stat) {
	overallStatus=stat;
	overallRatio = (float)stat.affected.count() / (float)stat.unaffected.count();
}

}

}

#endif
