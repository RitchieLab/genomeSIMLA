//
// C++ Interface: resultsrepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_REPORTINGRESULTSREPOSITORY_H
#define GENETICS_REPORTINGRESULTSREPOSITORY_H

#include "snpevaluation.h"
#include "ptestdistribution.h"
#include "modelstatistics.h"

namespace Genetics {

namespace Reporting {

using namespace Genetics::Evaluation;

/**
 * The statistics associated with a given cross validation slice
 */
typedef vector<FoldStatistics> SliceResults;		

/**
Accepts models and records the various statistics for final reporting. This report is designed to work within the cross validation framework taking the necessary averages

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ResultsRepository 
{
public:
    ResultsRepository(uint maxModelSize, uint foldCount);

    virtual ~ResultsRepository();

	/**
	 * Evaluation object stores each model it sees into the repository.
	 * @param order The size of the model being recorded
	 * @param slice The particular slice of the cross validation set
	 * @param stats The set of data associated with this entry
	 * @param label The label of the model being recorded
	 */
	virtual void Append(uint fold, FoldStatistics& stats);

	/**
	 * Assign that contains the distribution to be used during reporting. 
  	 * @note If no distribution is assigned or a NULL distribution is used, the p-value
	 * isn't reported. 
	 */
	void SetFitnessDistribution(PTestDistribution *dist);
	void SetOddsRatioDistribution( PTestDistribution *dist);
	void SetPEDistribution( PTestDistribution *dist );
	void SetLRDistribution( PTestDistribution *dist );
	/**
	 * Sorts the recorded objects and removes all but the top reportSize entries. 
	 */
	//virtual void Prune(uint reportSize);

	/**
	 * Stream report to the given stream.
	 * @param os The stream to be written to
	 * @param reportSize The number of models to report.
	 */
	virtual void GenerateReport(ostream* os, uint reportSize, uint start, uint stop);
	virtual void BakeResultModels(uint start, uint stop, uint reportSize, SnpRepository *repo, SnpEvaluation *eval);
	
	virtual void WriteHeader(ostream *os);
	/**
	 * Returns the model at position, position, for the given order. 
	 * @note This produces a model with a single "instance". So, callers should "Reduce" their copy. 
	 */
	//virtual FoldStatistics &GetModelStat(uint order, uint fold, uint position);

	/**
	 * Returns the statistics associated with position, position for all models of order size
	 */
	virtual ModelStatistics GetModelStatistics(uint order, uint position);


	void GenerateReport(ostream *os, uint start, uint stop, uint reportID, SnpRepository *repo, SnpEvaluation *eval);

	virtual void Sort();

protected:
	uint foldCount;							///<Used to record the number of folds in use
	uint maxModelSize;						///<The maximum number of SNPs in the models
	PTestDistribution *fitnessDist;			///<Distribution associated with the fitness value
	PTestDistribution *orDist;				///<Distribution associated with the odds ratio
	PTestDistribution *peDist;				
	PTestDistribution *lrDist;
	/**
	 * The results for each of the cross validation results. 
	 * @note Each entry in results represents a different cross validation "fold"
	 * @note result[order][fold]
	 */
	SliceResults *results;	
	bool isSorted;
	

	float GetAverageTesting(uint modelSize, uint position, float &oddsRatio, string& snpLabel);


};

inline
ResultsRepository::ResultsRepository(uint maxSize, uint folds) : foldCount(folds), maxModelSize(maxSize), fitnessDist(NULL), orDist(NULL), isSorted(false) { 
	results = new SliceResults[maxSize * foldCount];
	//results = new SliceResults[maxSize][foldCount];
}

inline
void ResultsRepository::SetPEDistribution( PTestDistribution *dist ) {
	peDist = dist;
}

inline
void ResultsRepository::SetLRDistribution (PTestDistribution *dist) {
	lrDist = dist;
}

inline
void ResultsRepository::SetOddsRatioDistribution( PTestDistribution *dist) {
	orDist = dist;
}

inline
void ResultsRepository::SetFitnessDistribution( PTestDistribution *dist) {
	fitnessDist = dist;
}

inline
void ResultsRepository::Sort() {
	if (!isSorted) {	
		for (uint i =0; i<maxModelSize; i++) {
			for (uint j=0; j<foldCount; j++) {
				uint idx = i * foldCount + j;
				sort(results[idx].begin(), results[idx].end());
			}
		}
		isSorted = true;
	}
}



inline
void ResultsRepository::Append(uint fold, FoldStatistics& stats) {
	isSorted = false;
	if (stats.order < maxModelSize && fold < foldCount) {
		uint idx = stats.order * foldCount + fold;
		assert(idx < maxModelSize * foldCount);
		results[idx].push_back(stats);
	}
	else {
		cout<<"Attempting to append a model that didn't fit inside the result repository\n";
		assert(0);
	}
}


}

}

#endif
