//
// C++ Implementation: orderreport
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "orderreport.h"
#include "evaluationmethod.h"

namespace MdrPDT {

namespace Evaluation {


OrderReport::OrderReport(int xvCount, int reportSize, int order) : 
		reportSize(reportSize), order(order) { 
	for (int i=0; i<xvCount; i++) {
		TReportTree tree(reportSize);
		reportFolds.push_back(tree);
	}
}

void OrderReport::AddModelResult(PdtModel& result) {
	int xvCount = reportFolds.size();
	for (int i=0; i<xvCount; i++) {
		float trainingT = result.GetTrainingT(i);
		TReportTree &report = reportFolds[i];
		report.Add(trainingT, result);
	}
}

/**
 * Skips forward from best to "idx"th position and returns the best model across all folds
 */
FoldStatistic OrderReport::GetBestModel(int idx) {
	int foldCount = reportFolds.size();
	float sumMOR = 0.0;

	FoldEvaluation foldEval;
	assert(idx < reportSize);
	for (int f=0; f<foldCount; f++) {
		TReportTreeNode *node = reportFolds[f].GetLast();

		int counter=0;
		//Skip ahead to the node of interest
		while (counter++ < idx)
			node = node->GetPrev();
		PdtModel report = node->GetData();
		if (foldCount == 1) 
			sumMOR += report.EvaluateMOR().GetRatio();
		else
			sumMOR += report.EvaluateMOR(f).GetRatio();

		foldEval.AddFold(report, f);
	}
	FoldStatistic fstat = foldEval.GetBestModel();
	fstat.SetAvgMOR(sumMOR/(float)foldCount);
	return fstat;		//foldEval.GetBestModel();
}

//For Todd's Biased MOR, we can find the bet, and then just return it
float OrderReport::GetAvgMOR(int idx) {
	int foldCount = reportFolds.size();
	float sumMOR = 0.0;
	FoldEvaluation foldEval;
	
	for (int f=0; f<foldCount; f++) {
		TReportTreeNode *node = reportFolds[f].GetLast();

		int counter=0;
		//Skip ahead to the node of interest
		while (counter++ < idx)
			node = node->GetPrev();
		PdtModel report = node->GetData();
		if (foldCount == 1) 
			sumMOR += report.EvaluateMOR().GetRatio();
		else
			sumMOR += report.EvaluateMOR(f).GetRatio();
	}
	return sumMOR / (float)foldCount;
}


void OrderReport::GenerateDetailedReport(int idx, std::ostream& os, EvaluationMethod *eval) {
	
	int foldCount = reportFolds.size();

	if (foldCount > 1) {
		for (int f=0; f<foldCount; f++) {
			if (reportFolds[f].GetLast()) {
				PdtModel model = reportFolds[f].GetLast()->GetData();
				eval->EvaluateModelVerbose(model, f, os);
			}
		}
	}
	else {
		FoldStatistic best = GetBestModel(idx);
		eval->EvaluateModelVerbose(best.model, os);		
	}
}

void OrderReport::GenerateSummaryReport(std::ostream& os, Distribution::PTestDistribution *distribution) {	
	std::vector<TReportTreeNode *> finalReport;
	int foldCount = reportFolds.size();

	for (int i=0; i<foldCount; i++)
		finalReport.push_back(reportFolds[i].GetLast());

	for (int m=0; m<reportSize; m++) {
		TReportTreeNode *node = finalReport[0];
		if (node) {
			FoldEvaluation eval;

			//Build up the Fold Evaluation 
			for (int i=0; i<foldCount; i++) {
				const PdtModel &report = finalReport[i]->GetData();
				eval.AddFold((PdtModel&)report, i);
				finalReport[i] = finalReport[i]->GetPrev();
			}
			if(m==0) {
				eval.ReportHeader(os, foldCount, distribution);
			}
			eval.ReportBestModel(os, m+1, foldCount, distribution);
		}
	}
}


}

}
