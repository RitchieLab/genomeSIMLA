//
// C++ Interface: foldstatistic
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTFOLDSTATISTIC_H
#define MDRPDTFOLDSTATISTIC_H
#include <vector>
#include "matchedoddsratio.h"
#include "ptestdistribution.h"
#include "pdtmodel.h"

namespace MdrPDT {

namespace Evaluation {

using namespace std;

struct floatCompare {
	int operator()(const float& l, const float& r) const {
		if (l<r) return -1;
		if (l>r) return 1;
		return 0;
	}
};


/**
	@author Eric Torstenson
*/
class FoldStatistic {
public:
	FoldStatistic();
	FoldStatistic(int totalFoldCount);
	int GetXVConsitency();
	
	void AddFold(PdtModel& model, int foldID);

	float AvgTrainingT() const;
	float AvgTestingT() const;

	//This requires a bit more work, so we'll only do it for the cases we care about
	float AvgMOR();
	void SetAvgMOR(float mor);
	float SumMOR();
	float AvgBiasedMOR();
	bool operator<(FoldStatistic& other);

	int GetCount();
		
	PdtModel model;
	vector<int> foldIDs;

	
//	float trainingT;
//	float testingT;

	int foldCount;					///<The number of cross validation folds there are
	float cvc;						///<Cross validation consistency
	float avgMOR;					///<The average MOR
	float sumMOR;					///<The sum of the mor
	float sumTrainingT;				///<The sum of the training T
	float sumTestingT;				///<The sum of the testing T
};
class FoldEvaluation {
public:
	void AddFold(PdtModel& report, int foldID);

	void ReportHeader(std::ostream& os, int foldCount, Distribution::PTestDistribution *dist);

	FoldStatistic GetBestModel();

	PdtModel ReportBestModel(std::ostream& os, int rank, int foldCount, Distribution::PTestDistribution *distribution);
	
	map<string, FoldStatistic> stats;
	vector<string> modelIDs;		///<The ID for each model associated with each fold

	static float reportThreshold;	///<Maximum pvalue for reporting
};



}

}

#endif
