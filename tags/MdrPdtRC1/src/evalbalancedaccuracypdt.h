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
#ifndef ESEEVALBALANCEDACCURACYPDT_H
#define ESEEVALBALANCEDACCURACYPDT_H
#include <math.h>

#include "genetics/snpevaluation.h"
#include "utility/utility.h"
#include <boost/timer.hpp>
#include "genetics/familyrepository.h"
namespace MDRPDT {

using namespace Utility;
using namespace Genetics;
using namespace boost;

struct MaxDPercentage {
	MaxDPercentage() : highestD(0.0), minD(0.0), percentage(0.65) {}
	MaxDPercentage(float perc) : highestD(0.0), minD(0.0), percentage(perc) {}
	
	/**
	 * @brief Evaluate currD to determine if it's within the threshold (updating the max D seen as needed)
	 * A D is passed in and it is compared with the minimum D. If the new value is below the minimum, we
	 * immediately return false. Otherwise, we keep up with the max observed so far and calculate the minimum
	 * according to our threshold. 
	 * @param currD the D being considered
	 */
	bool Evaluate(float currD) {
		if (currD < minD) 
			return false;
		if (percentage < 0.10)
			return true;
		if (currD > highestD) {
			highestD=currD;
			minD = highestD * percentage;
		}
		return true;
	}

	float highestD;					///<The highest D observed so far
	float minD;						///<The minimum D tolerable
	float percentage;				///<The threshold we want to consider to be valid
};

/**
 * @brief Used to refer to the type of status data is associated with a fold (for x-validation)
 */
typedef CaseControlStatus * FoldType;

/**
 * 	@brief This module is for evaluating a repository for unweighed balanced accuracy
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class EvalBalancedAccuracyPDT : public Genetics::Evaluation::SnpEvaluation {
public:
    class Results;
	
	/**
	 * @brief Construction
	 * @param comboStart The minimum size of models to be considered 
	 * @param comboEnd   The maximum size of models to be considered
	 * @param minDThresh The user assignable optimization value for D
	 * @param isPTest Indicate that this is the real run or a p-test (used to control verbosity)
	 * @param log The log where output will be written. NULL indicates that this won't generate output to a log
	 * @param testNumber The test number if this is a ptest
	 */
	EvalBalancedAccuracyPDT(uint comboStart, uint comboEnd, float minDThresh, bool isPTest, BasicLog *log = NULL, bool modelThresh=true, uint testNumber = 0);
	~EvalBalancedAccuracyPDT();

	/**
	 * @brief Adds a new fold to the set. 
	 * Each fold represents a set of weighted statuses to be used in identifying the top model
	 */
	void AppendStatus(FoldType fold, uint foldSize, FoldType pgStat, uint statSize);

	/**
	 * @brief Evaluate a model over the ptests and the status set. 
	 * This represents a single evaluation as opposed to a x-validation one	
	 */
	bool EvaluateModel(SnpAligned *snp, uint position);

	/**
	 * @brief Used for X-Validation, we will search through the folds for the best snp considering all model sizes.
	 */
	float GetTopScore(SnpAligned *&snp);
	/**
	 * @brief Used for X-Validation, we will search through the folds for best model for a given size
	 */
	float GetTopScore( uint modelSize, SnpAligned *&snp);
	
	/**
	 * @brief Perform final reporting
	 */
	void ReportResults(ostream *os);	

	/**
	 * @brief Perform Summary of the results
	 */
	void SummarizeResults(ostream *os);
		
	/**
	 * @brief Setup the data for final reporting. 
	 */
	void PostEvaluation();

	/**
	 * @brief Allow the application to force the description of specified models
	 */
	void EvaluateVerbose(SnpAligned *snp);
	
	/**
	 * @brief Reports the Matched Odds Ratio for a given model (and High Risk Cells)
	 * @param snp The snp being considered
	 * @param status The array of parental statuses being considered
	 * @param highRiskCells The array of genotype indices that are deemed HR
	 * @param hrCellCount The number of entries in the HR cell array
	 */
	void FindMatchedOddsRatio(SnpAligned *snp, FoldType status, uint highRiskCells[], uint hrCellCount);

	static bool Verbose;					///<Used to toggle the production of all T values during evaluation

	/**
	 * @brief Evaluate a given snp for a specified status 
	 */
	float Evaluate(SnpAligned *, FoldType fold);

protected:

	/**
	 * @brief Perform verbose evaluation 
	 */
	float EvaluateVerbose(SnpAligned *, FoldType fold);

	/**
	 * @brief Array of status arrays
	 */	
	vector <FoldType> testingSet;				///<Set of vectors associated with the current fold
	vector<FoldType> sibshipStatus;				///<Vector containing the siships associated with the current data set

	vector <Results *> results;					///<Results are stored for each of the folds
	bool isPTest;
	uint comboEnd;								///<Max loci of a model
	uint comboStart;							///<Min loci count of a model
	uint testNumber;
	MaxDPercentage *qdMax;						///<The array of quick D maxes observed

	BasicLog *log;								///<Destination for report log
	uint statusCount;							///<The size of the testingSet array
	uint sibshipStatusCount;					///<The size of the sibshipStatus array

public:
/**
 * @brief Results associated with a single BA run
 */
struct Results {
	Results(uint comboStart, uint comboEnd) : comboEnd(comboEnd), comboStart(comboStart), topModel(comboEnd, NULL), highTStatistic(comboEnd, 0.0) {}
	~Results() { 
		for (uint i=comboStart-1; i<comboEnd; i++) 
			if (topModel[i])
				topModel[i]->ReduceInstanceCount(); 
		}
	
	/**
	 * @brief Perform evaluation for a specific model
	 */		
	void Evaluate(SnpAligned *snp, float score);
	uint comboEnd;
	uint comboStart;
	vector<SnpAligned *>topModel;			///<The top model observed so far
	vector<float>highTStatistic;			///<The winning t-statistic
};
};



inline
void EvalBalancedAccuracyPDT::Results::Evaluate( SnpAligned *snp, float score) {
	uint modelSize=snp->GetLabelCount()-1;
	if (modelSize < highTStatistic.size() && score > highTStatistic[modelSize]) {
		highTStatistic[modelSize] = score;
		SnpAligned *curTop = topModel[modelSize];
		if (curTop)
			curTop->ReduceInstanceCount();
		snp->IncrementInstanceCount();
		topModel[modelSize] = snp;
	}
}


inline
float EvalBalancedAccuracyPDT::GetTopScore( uint modelSize, SnpAligned *&snp) {
	uint count=results.size();
	float bestResult =0.0;
	for (uint i=0; i<count; i++) {
		if (bestResult < results[i]->highTStatistic[modelSize]) {
			bestResult=results[i]->highTStatistic[modelSize];
			snp=results[i]->topModel[modelSize];
		}
	}
	return bestResult;
}

// TODO I need to set this up to work with x-validation. 
// Currently, we are looking for the over winner!
inline
float EvalBalancedAccuracyPDT::GetTopScore(SnpAligned *&snp) {
	uint count=results.size();
	float bestResult =0.0;
	for (uint i=0; i<count; i++) {
		for (uint p=comboStart-1; p<comboEnd; p++) { 
			if (bestResult < results[i]->highTStatistic[p]) {
				bestResult=results[i]->highTStatistic[p];
				snp=results[i]->topModel[p];
			}
		}
	}
	return bestResult;
}

inline
EvalBalancedAccuracyPDT::EvalBalancedAccuracyPDT(uint comboStart, uint comboEnd, float minDThresh, bool isPTest, BasicLog *log /*= NULL*/, bool modelThresh, uint testNumber/*=0*/) : SnpEvaluation(modelThresh), isPTest(isPTest), comboEnd(comboEnd), comboStart(comboStart), testNumber(testNumber), log(log), statusCount(0), sibshipStatusCount(0)  {
	qdMax=new MaxDPercentage[comboEnd];
	for (uint i=0; i<comboEnd; i++)
		qdMax[i].percentage=minDThresh;
}

inline
EvalBalancedAccuracyPDT::~EvalBalancedAccuracyPDT() {
	uint count=results.size();
	for (uint i=0; i<count; i++) 
		delete results[i];
	results.clear();
	if (qdMax)
		delete[] qdMax;
}


inline
void EvalBalancedAccuracyPDT::AppendStatus( FoldType fold, uint foldSize, FoldType pgFold, uint pgFoldSize) {
	if (statusCount != 0) 
		assert(statusCount == foldSize);
	statusCount = foldSize;
	testingSet.push_back(fold);

	if (pgFold) {
		if (sibshipStatusCount != 0) {
			assert(sibshipStatusCount == foldSize);
		}
		sibshipStatusCount = foldSize;
		sibshipStatus.push_back(fold);

	}

	//testingSet=fold;
	Results *newResults = new Results(comboStart, comboEnd);
	results.push_back(newResults);
}


inline
bool EvalBalancedAccuracyPDT::EvaluateModel(SnpAligned *snp, uint position) {
	uint count=results.size();
	float score=0.0;


	//Iterate over each of the status objects insde the vector
	for (uint i=0; i<count; i++) {
		FoldType ccSet = testingSet[i];
		score=Evaluate(snp, ccSet);
		results[i]->Evaluate(snp, score);
	}
	return score > 0.0;
}

inline
void EvalBalancedAccuracyPDT::EvaluateVerbose(SnpAligned *snp) {
	uint count = testingSet.size();
	//Iterate over each of the status objects insde the vector
	for (uint i=0; i<count; i++) {
		FoldType ccSet = testingSet[i];
		EvaluateVerbose(snp, ccSet);
	}
}




}

#endif
