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
#include <math.h>

namespace Genetics {

namespace Evaluation {

using namespace Genetics;

/**
 * This is required to build up the rest of the statistics
 */
void ModelStatistics::GatherTestStatistics(Genetics::SnpRepository *repo, SnpEvaluation *eval, ostream *os) {
	stringstream xvc;
	map<string, TestHolder> frequencies;
	uint count 		= folds.size();
	sumTestingStatistic = 0.0;
	sumPredError 	= 0.0;
	sumClassError 	= 0.0;
	sumMOR			= 0.0;

	if (count < 1) {
		cout<<"Fold count for statistic fold list is 0. \n";
		return;
	}


/*	cout<<"     "<<setw(5)<<"Ord."<<setw(5)<<"fold"<<setw(12)<<"Label"
		<<setw(12)<<"Training"<<setw(12)<<"Cls."
		<<setw(12)<<"Testing"<<setw(12)<<"Prd."
		<<setw(12)<<"MOR"<<"\n";
*/
	for (uint i=0; i<count; i++) {
		FoldStatistics &fold = folds[i];
		
		SnpAligned *snp = repo->GetSnp(fold.label.c_str());
		fold = eval->GetFoldStatistics(i, snp);
		frequencies[fold.label].AddTest(fold);
	
		if (isnormal(fold.testingFitness)) {
			sumTestingStatistic+=fold.testingFitness;
			testingCount++;
		}

		if (isnormal(fold.predError)) {
			sumPredError+=fold.predError;
			predCount++;
		}

		if (isnormal(fold.clsError)) {
			sumClassError+=fold.clsError;
			clsCount++;
		}

		if (isnormal(fold.oddsRatio)) {
			sumMOR+=fold.oddsRatio;
			orCount++;
		}

		if (isnormal(fold.trainingFitness)) {
			sumTrainingStatistic+=fold.trainingFitness;
			trainingCount++;
		}

/*		cout<<"fold "<<i<<": ";
		cout<<setw(5)<<fold.order<<setw(5)<<i<<setw(12)<<fold.label;
		cout<<setw(12)<<setprecision(5)<<fold.trainingFitness;
		cout<<setw(12)<<setprecision(2)<<fold.clsError;
		cout<<setw(12)<<setprecision(5)<<fold.testingFitness;
		cout<<setw(12)<<setprecision(2)<<fold.predError;
		cout<<setw(12)<<setprecision(5)<<fold.oddsRatio<<"\n";
*/
		if (snp)
			snp->ReduceInstanceCount();
		if (os) {
			*os<<setw(16)<<fold.order<<setw(16)<<i<<setw(12)<<fold.label;
			*os<<setw(12)<<setprecision(5)<<fold.trainingFitness;
			*os<<setw(12)<<setprecision(2)<<fold.clsError;
			*os<<setw(12)<<setprecision(5)<<fold.testingFitness;
			*os<<setw(12)<<setprecision(2)<<fold.predError;
			*os<<setw(12)<<setprecision(5)<<fold.oddsRatio<<"\n";
		}
	}

	vector<TestHolder> winners;					///<We'll stash each of the winners to find the real winner
	map<string, TestHolder>::iterator itr = frequencies.begin();
	map<string, TestHolder>::iterator end = frequencies.end();
	
//	int ack=0;
	for (; itr!=end; itr++) {
/*		cout<<"\t"<<setw(12)<<itr->second.snpLabel;
		cout<<setw(12)<<setprecision(5)<<itr->second.fitness;
		cout<<setw(12)<<setprecision(2)<<itr->second.freq;
		cout<<setw(12)<<setprecision(5)<<itr->second.oddsRatio<<"\n";
*/
		winners.push_back(itr->second);
	}
	
	sort(winners.begin(), winners.end());
	int topModel = winners.size() - 1;

/*Debugging tool
	for (uint ack=0; ack<winners.size(); ack++) {
		if (ack == topModel)
			cout<<"*";
		cout<<"\t"<<setw(12)<<winners[ack].snpLabel;
		cout<<setw(12)<<setprecision(5)<<winners[ack].fitness;
		cout<<setw(12)<<setprecision(2)<<winners[ack].freq;
		cout<<setw(12)<<setprecision(5)<<winners[ack].oddsRatio<<"\n";
	}
*/

	TestHolder &winner = winners[topModel];

	label=winner.snpLabel;	
	xvConst = winner.freq;




	//oddsRatio = winners[winners.size() - 1].oddsRatio;
	xvc<<xvConst<<"/"<<count;
	xvConsistency = xvc.str();

}

void ModelStatistics::Reset() {
	foldCount 		= 0;
	sumTestingStatistic = 0.0;
	sumTrainingStatistic = 0.0;
	sumPredError 	= 0.0;
	sumMOR			= 0.0;
	sumClassError 	= 0.0;


	testingCount	= 0;
	predCount 		= 0;
	clsCount 		= 0;
	orCount			= 0;
	trainingCount	= 0;

	//oddsRatio 		= 0.0;
	winnerFrequency = 0;
	folds.clear();
}

void ModelStatistics::Append(uint fold, FoldStatistics& stats) {
	sumTrainingStatistic+=stats.trainingFitness;
	folds.push_back(stats);
	foldCount++;
}
	


}

}
