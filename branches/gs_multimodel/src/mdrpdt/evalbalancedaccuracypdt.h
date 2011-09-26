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
#include "pdtfold.h"
#include "genetics/resultsrepository.h"
#include "genetics/snprepository.h"
#include "genetics/modelstatistics.h"

namespace MDR {

using namespace Utility;
using namespace Genetics;
using namespace boost;
using namespace Genetics::Evaluation;
using namespace Genetics::Reporting;

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
	 * @brief Configure the folds (array of size foldCount) with status of pgFold
	 */
	void SetFolds( PdtFold *fold, uint foldCount, FoldType pgFold, uint pgFoldSize);
	
	/**
	 * @brief Returns a bitset representing high risk cells (of size snp->CountGenotypes())
	 * @param snp The model to be evaluated
	 * @param individuals Is the set of individuals associated with the query. 
	 * @note individuals can be all or a partial segment (such as a fold)
	 */
	BitSetType GetHighRiskCells(SnpAligned *snp, CaseControlStatus &individuals);

	/**
	 * @brief Returns a bitset representing high risk cells (of size snp->CountGenotypeS())
	 * @param snp The model to be evaluated
	 * @note This will return the high risk for all individuals with valid data
	 */
	BitSetType GetHighRiskCells( SnpAligned * snp);
	
	void ReportResults(ostream *os, uint reportSize, SnpRepository *repo);
	void BakeResultModels(uint reportSize, SnpRepository *repo);
	/**
	 * @brief Evaluate a model over the ptests and the status set. 
	 * This represents a single evaluation as opposed to a x-validation one	
	 */
	bool EvaluateModel(SnpAligned *snp, uint position);

	/**
	 * @brief Used for X-Validation, we will search through the folds for the best snp considering all model sizes.
	 */
	ModelStatistics GetModelStats(SnpAligned *&snp);
	/**
	 * @brief Used for X-Validation, we will search through the folds for best model for a given size
	 */
	ModelStatistics GetModelStats(uint modelSize, SnpAligned *&snp);
	FoldStatistics GetFoldStatistics( uint modelSize, SnpAligned *snp);
	
	float GetAverageTesting(uint modelSize, float &oddsRatio, SnpAligned*& snp);
	float GetAverageTesting(SnpAligned *& snp);

	/**
	 * @brief Perform final reporting
	 */
	void ReportResults(ostream *os);	

	//void ReportResults(ostream *os, uint reportSize);

	/**
	 * @brief Perform Summary of the results
	 */
	void SummarizeResults(ostream *os);
		
	/**
	 * @brief Setup the data for final reporting. 
	 */
	void PostEvaluation();

	/**
	 * @brief Write the findings to the report (this is detailed findings, not the distribution findings)
	 */
	void DisplayAnalyses(bool writeDetails);

	/**
	 * @brief Allow the application to force the description of specified models
	 */
	ModelStatistics EvaluateVerbose(SnpAligned *snp);
	FoldStatistics EvaluateVerbose(uint fold, SnpAligned *snp, stringstream& details);
//	FoldStatistics EvaluateFold(uint fold, SnpAligned *snp);

	/**
	 * @brief Reports the Matched Odds Ratio for a given model (and High Risk Cells)
	 * @param snp The snp being considered
	 * @param status The array of parental statuses being considered
	 * @param highRiskCells The array of genotype indices that are deemed HR
	 * @param hrCellCount The number of entries in the HR cell array
	 */
	float FindMatchedOddsRatio(SnpAligned *snp, FoldType status);
	float FindMatchedOddsRatio(SnpAligned *snp, FoldType status, CaseControlStatus &mask);
	static bool Verbose;					///<Used to toggle the production of all T values during evaluation

	void SetFitnessDist(PTestDistribution *dist);
	void SetOddsRatioDist(PTestDistribution *dist);
	void SetPEDist(PTestDistribution *dist);
protected:

	/**
	 * @brief Array of status arrays
	 */	
	PdtFold *folds;								///<Array of cross validation slices
	uint foldCount;								///<The number of cross validation slices
	//vector <FoldType> testingSet;				///<Set of vectors associated with the current fold
	vector<FoldType> sibshipStatus;				///<Vector containing the siships associated with the current data set

	Results *results;							///<Results are stored for each of the folds
	bool isPTest;
	uint comboEnd;								///<Max loci of a model
	uint comboStart;							///<Min loci count of a model
	uint testNumber;
	MaxDPercentage *qdMax;						///<The array of quick D maxes observed

	BasicLog *log;								///<Destination for report log
	uint statusCount;							///<The size of the testingSet array
	uint sibshipStatusCount;					///<The size of the sibshipStatus array

	ResultsRepository *overallResults;		 	///<results used in overall report

public:



/**
 * @brief Results associated with a single BA run
 */
struct Results {

	Results() : comboEnd(0), comboStart(0) {}
	Results(uint comboStart, uint comboEnd) : comboEnd(comboEnd), comboStart(comboStart), topModel(comboEnd, NULL), topStatistics(comboEnd) {} 
	//, highTStatistic(comboEnd, 0.0), classError(comboEnd, 0.0), matchedOddsRatio(comboEnd, 0.0) {}
	
	void SetSize(uint comboStart, uint comboEnd) {
		this->comboStart = comboStart;
		this->comboEnd = comboEnd;
		topModel = vector<SnpAligned *>(comboEnd, NULL);
		topStatistics = vector<FoldStatistics>(comboEnd);
		//highTStatistic = vector<float>(comboEnd, 0.0);
		//classError     = vector<float>(comboEnd, 0.0);
		//matchedOddsRatio = vector<float>(comboEnd, 0.0);
		
	}
	
	~Results() { 
		for (uint i=comboStart-1; i<comboEnd; i++) 
			if (topModel[i])
				topModel[i]->ReduceInstanceCount(); 
		}
	
	/**
	 * @brief Perform evaluation for a specific model
	 */		
	//void Evaluate(SnpAligned *snp, PdtFold::AccuracyEval &score);
	void Evaluate(SnpAligned *snp, FoldStatistics &score);
	uint comboEnd;
	uint comboStart;
	vector<SnpAligned *>topModel;			///<The top model observed so far
	vector<FoldStatistics> topStatistics;				///<The scores associated with the various pieces
	//vector<float>highTStatistic;			///<The winning t-statistic
	//vector<float>classError;				///<The classification error for the winner
	//vector<float>matchedOddsRatio;			///<The Matched odds ratio for the winner
};
};


inline
void EvalBalancedAccuracyPDT::SetPEDist(PTestDistribution *dist) {
	peDist = dist;
	if (overallResults)
		overallResults->SetPEDistribution(dist);
}

inline
void EvalBalancedAccuracyPDT::SetFitnessDist(PTestDistribution *dist) {
	fitnessDist=dist;
	if (overallResults)
		overallResults->SetFitnessDistribution(dist);
}


inline
void EvalBalancedAccuracyPDT::SetOddsRatioDist(PTestDistribution *dist) {
	orDist = dist;
	if (overallResults)
		overallResults->SetOddsRatioDistribution(dist);
}


inline
//void EvalBalancedAccuracyPDT::Results::Evaluate( SnpAligned *snp, PdtFold::AccuracyEval &score) {
void EvalBalancedAccuracyPDT::Results::Evaluate( SnpAligned *snp, FoldStatistics &score) {
	uint modelSize=snp->GetLabelCount()-1;
	if (modelSize < topStatistics.size() && score.trainingFitness > topStatistics[modelSize].trainingFitness) {
		topStatistics[modelSize] = score;
		SnpAligned *curTop = topModel[modelSize];
		if (curTop)
			curTop->ReduceInstanceCount();
		snp->IncrementInstanceCount();
		topModel[modelSize] = snp;
		
	}
}


struct TestHolder {
	int freq;
	float fitness;
	float oddsRatio;
	SnpAligned *snp;
	
	TestHolder() : freq(0), fitness(0.0), snp(NULL) { }
	
	void AddTest(SnpAligned *snp, FoldStatistics &stats) {
		if (this->snp) 					{ assert(snp == this->snp); 	}

		freq++;
		this->snp = snp;
		this->fitness+=stats.testingFitness;
		this->oddsRatio=stats.oddsRatio;
	}

	bool operator<(const TestHolder& other) const  { 	return (freq <= other.freq) || (freq == other.freq && AverageFitness() < other.AverageFitness()); 	}

	bool operator==(TestHolder& other) { return freq == other.freq || AverageFitness() == other.AverageFitness(); }
	float AverageFitness() const {    			return fitness/(float)freq; 	}		

};


inline
float EvalBalancedAccuracyPDT::GetAverageTesting(uint modelSize, float &oddsRatio, SnpAligned*& snp) {
	map<string, TestHolder> frequencies;
	for (uint i=0; i<foldCount; i++) {
		SnpAligned *snp = results[i].topModel[modelSize];
		frequencies[snp->GetLabel()].AddTest(snp, results[i].topStatistics[modelSize]);
	}

	//TestHolder winners[frequencies.size()];
	vector<TestHolder> winners;
	map<string, TestHolder>::iterator itr = frequencies.begin();
	map<string, TestHolder>::iterator end = frequencies.end();
	
	for (; itr!=end; itr++) 
		winners.push_back(itr->second);
		//winners[i++] = itr->second;

	sort(winners.begin(), winners.end());
	float value = winners[winners.size() - 1].AverageFitness();
	snp=winners[winners.size() - 1].snp;	
	oddsRatio = winners[winners.size() - 1].oddsRatio;
	return value;

}

inline
float EvalBalancedAccuracyPDT::GetAverageTesting(SnpAligned*& snp) {
	float value = 0;
	cout<<"We don't know how to pick the best model over different model sizes\n";
	assert(0);
	return value;

}


inline
ModelStatistics EvalBalancedAccuracyPDT::GetModelStats( uint modelSize, SnpAligned *&snp) {
	ModelStatistics stats;
	FoldStatistics bestResult("");
	for (uint i=0; i<foldCount; i++) {
		Results &r = results[i];
		FoldStatistics& foldData = r.topStatistics[modelSize];
		stats.Append(i, results[i].topStatistics[modelSize]);
	}
	return stats;
}

// TODO I need to set this up to work with x-validation. 
// Currently, we are looking for the over winner!
inline
ModelStatistics EvalBalancedAccuracyPDT::GetModelStats(SnpAligned *&snp) {
	FoldStatistics bestResult("");
	ModelStatistics finalResults;
		for (uint p=comboStart-1; p<comboEnd; p++) { 
			ModelStatistics curStats = GetModelStats(snp);
			assert(0);
/*			if (bestResult < results[i].topStatistics[p]) {
				bestResult=results[i].topStatistics[p];
				snp=results[i].topModel[p];
				bestResult.label = snp->GetLabel();
			}
		} */
	}
	return finalResults;
}

inline
EvalBalancedAccuracyPDT::EvalBalancedAccuracyPDT(uint comboStart, uint comboEnd, float minDThresh, bool isPTest, BasicLog *log /*= NULL*/, bool modelThresh, uint testNumber/*=0*/) : SnpEvaluation(modelThresh), folds(NULL), isPTest(isPTest), comboEnd(comboEnd), comboStart(comboStart), testNumber(testNumber), log(log), statusCount(0), sibshipStatusCount(0), overallResults(NULL)  {
	qdMax=new MaxDPercentage[comboEnd];
	for (uint i=0; i<comboEnd; i++)
		qdMax[i].percentage=minDThresh;
}

inline
EvalBalancedAccuracyPDT::~EvalBalancedAccuracyPDT() {
	if (results)
		delete[] results;

	if (qdMax)
		delete[] qdMax;

	if (overallResults)
		delete overallResults;

	if (folds)
		delete[] folds;
}


inline
void EvalBalancedAccuracyPDT::SetFolds( PdtFold *fold, uint foldCount, FoldType pgFold, uint pgFoldSize) {
	if (folds)
		delete[] folds;

	folds = new PdtFold[foldCount];
	//Deep copy the folds
	for (uint i=0; i<foldCount; i++) 
		folds[i] = fold[i];

//	folds=fold;
	
	this->foldCount = foldCount;
	results = new Results[foldCount];

	if (overallResults)
		delete overallResults;
	overallResults = new ResultsRepository(comboEnd, foldCount);

	for (uint i=0; i<foldCount; i++) 
		results[i].SetSize(comboStart, comboEnd);

	//This looks suspicious to me
	if (pgFold) {
		if (sibshipStatusCount != 0) {
			assert(sibshipStatusCount == pgFoldSize);
		}
		sibshipStatusCount = pgFoldSize;
		sibshipStatus.push_back(pgFold);
	}

}







}

#endif
