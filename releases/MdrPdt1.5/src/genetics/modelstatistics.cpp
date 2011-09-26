//
// C++ Implementation: modelstatistics
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genetics/snprepository.h"
#include "modelstatistics.h"
#include "genetics/snpevaluation.h"
#include <iomanip>

namespace Genetics {

namespace Evaluation {

using namespace Genetics;

/**
 * This is required to build up the rest of the statistics
 */
void ModelStatistics::GatherTestStatistics(Genetics::SnpRepository *repo, SnpEvaluation *eval, ostream *os) {
	stringstream xvc;
	map<string, TestHolder> frequencies;
	uint count = folds.size();
	sumTestingStatistic = 0.0;
	sumPredError = 0.0;
	sumClassError = 0.0;
	for (uint i=0; i<count; i++) {
		FoldStatistics &fold = folds[i];
		
		SnpAligned *snp = repo->GetSnp(fold.label.c_str());
		fold = eval->GetFoldStatistics(i, snp);
		frequencies[fold.label].AddTest(fold);
		sumTestingStatistic+=fold.testingFitness;
		sumPredError+=fold.predError;
		sumClassError+=fold.clsError;
		snp->ReduceInstanceCount();
		if (os) {
			*os<<setw(16)<<fold.order<<setw(16)<<i<<setw(12)<<fold.label;
			*os<<setw(12)<<setprecision(3)<<fold.trainingFitness;
			*os<<setw(12)<<setprecision(2)<<fold.clsError;
			*os<<setw(12)<<setprecision(3)<<fold.testingFitness;
			*os<<setw(12)<<setprecision(2)<<fold.predError;
			*os<<setw(12)<<setprecision(3)<<fold.oddsRatio<<"\n";
		}
	}

	vector<TestHolder> winners;					///<We'll stash each of the winners to find the real winner
	map<string, TestHolder>::iterator itr = frequencies.begin();
	map<string, TestHolder>::iterator end = frequencies.end();
	
	for (; itr!=end; itr++) 
		winners.push_back(itr->second);
		//winners[i++] = itr->second;

	sort(winners.begin(), winners.end());
	label=winners[winners.size() - 1].snpLabel;	
	oddsRatio = winners[winners.size() - 1].oddsRatio;
	xvc<<winners[winners.size() - 1].freq<<"/"<<count;
	xvConsistency = xvc.str();
}

void ModelStatistics::Reset() {
	foldCount = 0;
	sumTestingStatistic = 0.0;
	sumTrainingStatistic = 0.0;
	sumPredError = 0.0;
	sumClassError = 0.0;
	folds.clear();
}

void ModelStatistics::Append(uint fold, FoldStatistics& stats) {
	sumTrainingStatistic+=stats.trainingFitness;
	folds.push_back(stats);
	foldCount++;
}
	


}

}
