//
// C++ Interface: evalbalancedaccuracy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEEVALBALANCEDACCURACY_H
#define ESEEVALBALANCEDACCURACY_H

#include "genetics/snpevaluation.h"
#include "utility/utility.h"
#include "genetics/modelstatistics.h"


namespace MDR {

using namespace Utility;
using namespace Genetics;



/**
 * 	@brief This module is for evaluating a repository for unweighed balanced accuracy
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class EvalBalancedAccuracy : public Genetics::Evaluation::SnpEvaluation {
public:
    class Results;
	
	/**
	 * @brief Construction
	 * @param eval Pointer to the actual evaluation method.
	 */
	EvalBalancedAccuracy();
	EvalBalancedAccuracy(uint loci) : highBA(0) {}

	~EvalBalancedAccuracy();

	/**
	 * @brief Adds a new fold to the set. 
	 * Each fold represents a set of weighted statuses to be used in identifying the top model
	 */
	void AppendStatus(CaseControlStatus &fold);

	/**
	 * @brief Evaluate a model over the ptests and the status set. 
	 * This represents a single evaluation as opposed to a x-validation one	
	 */
	bool EvaluateModel(SnpAligned *snp, uint position);

	Genetics::Evaluation::FoldStatistics GetTopScore(SnpAligned *&snp);
	float GetAverageTesting(SnpAligned *& snp) { assert(0); return 0.0;}

	/**
	 * @brief Used for X-Validation, we will search through the folds for the best snp considering all model sizes.
	 */
	Genetics::Evaluation::ModelStatistics GetModelStats(SnpAligned *&snp)	{ assert(0); return Genetics::Evaluation::ModelStatistics(); }
	/**
	 * @brief Used for X-Validation, we will search through the folds for best model for a given size
	 */
	Genetics::Evaluation::ModelStatistics GetModelStats(uint modelSize, SnpAligned *&snp)	{ assert(0); return Genetics::Evaluation::ModelStatistics();}
	Genetics::Evaluation::FoldStatistics GetFoldStatistics( uint modelSize, SnpAligned *snp)	{ assert(0); return Genetics::Evaluation::FoldStatistics();}

	void ReportResults(ostream *os);
protected:
	/**
	 * @brief Evaluate a given snp for a specified status 
	 */
	float Evaluate( SnpAligned *snp, const CaseControlStatus &status);

	/**
	 * @brief The array of tests
	 */
	vector <CaseControlStatus> testingSet;

	/**
	 * @brief The array of result objects
	 */
	vector <Results *> results;
	
public:
	//Quick test variable
	float highBA;
/**
 * @brief Results associated with a single BA run
 */
struct Results {
	Results() : topModel(NULL), statistics(NULL) {}
	~Results() {
		topModel->ReduceInstanceCount();
	}
	/**
	 * @brief Perform evaluation of a model to determine if it should be recorded 
	 */
	void Evaluate(SnpAligned *snp, Genetics::Evaluation::FoldStatistics stat) {
		if (statistics < stat) {
			statistics = stat;
			if (topModel)
				topModel->ReduceInstanceCount();
			topModel = snp;
			topModel->IncrementInstanceCount();
		}
	}
			
	SnpAligned *topModel;						///<The top model associated with this particular Result object
	Genetics::Evaluation::FoldStatistics statistics;								///<The score associated with the top model
};
};

// TODO I need to set this up to work with x-validation
inline
Genetics::Evaluation::FoldStatistics EvalBalancedAccuracy::GetTopScore(SnpAligned *&snp) {
	uint count=results.size();
	Genetics::Evaluation::FoldStatistics bestResult(NULL);
	for (uint i=0; i<count; i++) { 
		if (bestResult < results[i]->statistics) {
			bestResult=results[i]->statistics;
			snp=results[i]->topModel;
			bestResult.label=snp->GetLabel();
		}
	}
	return bestResult;
}
inline
EvalBalancedAccuracy::EvalBalancedAccuracy()  {
	
}

inline
EvalBalancedAccuracy::~EvalBalancedAccuracy() {
	uint count=results.size();
	for (uint i=0; i<count; i++) 
		delete results[i];
	results.clear();
}

inline
void EvalBalancedAccuracy::AppendStatus( CaseControlStatus &fold) {
	testingSet.push_back(fold);
	Results *newResults = new Results();
	results.push_back(newResults);
}


inline
bool EvalBalancedAccuracy::EvaluateModel(SnpAligned *snp, uint position) {
	uint count=results.size();

	for (uint i=0; i<count; i++) {
		Genetics::Evaluation::FoldStatistics score(snp->GetLabel());
		score.trainingFitness = Evaluate(snp, testingSet[i]);
		results[i]->Evaluate(snp, score);
	}
	return true;
}

inline
float EvalBalancedAccuracy::Evaluate( SnpAligned *snp, const CaseControlStatus &status) {
	float affCount = 0.0;										///<The number of bits masked to affected
	float unaffCount = 0.0;
	float ba=0.0;

	BitSetType curGenotype=snp->GetGenotype(0);
	uint gtCount = snp->CountGenotypes();


	float caseCount = (float)(curGenotype & status.affected).count(); 
	float controlCount = (float)(curGenotype & status.unaffected).count();
	float sensitivity = 0.0;
	float specificity = 0.0;
	float threshold = caseCount / controlCount;

	if (caseCount == 0 && controlCount == 0)
		return 0.0; 

	for (uint i=1; i<gtCount; i++) {
		curGenotype=snp->GetGenotype(i);
		affCount = (float)(curGenotype & status.affected).count();
		unaffCount = (float)(curGenotype & status.unaffected).count();

		//Mark as a High Risk
		if (unaffCount == 0 || affCount / unaffCount >= threshold) 
			sensitivity+=affCount;
		else
			specificity+=unaffCount;
	}		

	sensitivity/=caseCount;
	specificity/=controlCount;
	ba=((sensitivity + specificity))*(50.0);
	//cout<<snp->GetLabel()<<" "<<ba<<"\n";
	return ba;
}
}

#endif
