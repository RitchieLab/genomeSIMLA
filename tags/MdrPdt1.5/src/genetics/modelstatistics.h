//
// C++ Interface: modelstatistics
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSMODELSTATISTICS_H
#define GENETICSMODELSTATISTICS_H

#include "genetics/snpaligned.h"
#include <vector>
#include "snprepository.h"

namespace Genetics {

namespace Evaluation {

using namespace Genetics;
class SnpEvaluation;


using namespace std;



/**
 * @brief Cache for the various results
 */
struct FoldStatistics {
	float trainingFitness;		///<This is the metric associated with training
	float testingFitness;		///<This is the metric associated with testing

	float clsError;				///<This is the error (i.e. classification error)
	float predError;			///<This is the prediction error
	float oddsRatio;			///<The odds ratio
	string label;				///<The label associated with the model that produced these results
	uint order;					///<The order associated with the model that produced these results
	
	FoldStatistics(const char *label, uint order, float training, float testing, float clsError, float predError, float oddsRatio) : 
		trainingFitness(training), 
		testingFitness(testing), 
		clsError(clsError), 
		predError(predError), 
		oddsRatio(oddsRatio), 
		label(label), 
		order(order) {}

	FoldStatistics(const char *label = "") : trainingFitness(0.0), testingFitness(0.0), clsError(0.0), predError(0.0), oddsRatio(0.0), label(label), order(0) {}

	bool operator<(const FoldStatistics& other) const {
		return trainingFitness>other.trainingFitness;
	}
};


struct TestHolder {
	int freq;
	float fitness;
	float oddsRatio;
	string snpLabel;
	
	TestHolder() : freq(0), fitness(0.0), snpLabel("") { }
	
	void AddTest(FoldStatistics &stats) {
		assert(snpLabel == "" || strcmp(stats.label.c_str(), snpLabel.c_str()) == 0);  	
		freq++;
		snpLabel = stats.label;
		fitness+=stats.testingFitness;
		oddsRatio=stats.oddsRatio;
	}

	bool operator<(const TestHolder& other) const  { 	return freq < other.freq || AverageFitness() < other.AverageFitness(); 	}
	bool operator==(TestHolder& other) { return freq == other.freq || AverageFitness() == other.AverageFitness(); }
	float AverageFitness() const {    			return fitness/(float)freq; 	}		

};

/**
Holds the statistics for an entire model

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ModelStatistics {
public:
    ModelStatistics() : foldCount(0), sumTestingStatistic(0.0), sumTrainingStatistic(0.0), sumPredError(0.0), sumClassError(0.0), oddsRatio(0.0) {}
    ~ModelStatistics() {}

	/**
	 * This is required to build up the rest of the statistics
	 */
	void GatherTestStatistics(Genetics::SnpRepository *repo, SnpEvaluation *eval, ostream *os);

	void Reset();

	void Append(uint fold, FoldStatistics& stats);
	
	float GetAvgTraining() {
		return sumTrainingStatistic / (float)foldCount;
	}

	float GetAvgTesting() {
		return sumTestingStatistic / (float)foldCount;
	}
	
	float GetAvgPredictionError() {
		return sumPredError / (float)foldCount;
	}

	float GetAvgClassificationError() {
		return sumClassError / (float)foldCount;
	}
	
	float GetOddsRatio() {
		return oddsRatio;
	}

	string GetLabel() {
		return label;
	}

	vector<FoldStatistics> folds;
	uint foldCount;
	float sumTestingStatistic;
	float sumTrainingStatistic;
	float sumPredError;
	float sumClassError;
	float oddsRatio;
	string xvConsistency;
	string label;

};

}

}

#endif
