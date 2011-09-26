//
// C++ Implementation: tfinalreport
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tfinalreport.h"

namespace MdrPDT {

namespace Evaluation {
bool TFinalReport::VerboseAnalysis = false;


TFinalReport::TFinalReport(int maxModelSize, int xvCount, int reportSize) : reportSize(reportSize) {
	for (int i=0; i<maxModelSize; i++) {
		OrderReport report(xvCount, reportSize, i);
		reports.push_back(report);
	}
}

TFinalReport::~TFinalReport() { }

void TFinalReport::AddModelResult(PdtModel& result) {
	int order = result.ModelSize() - 1;
	reports[order].AddModelResult(result);
	if (VerboseAnalysis)
		cout<<result.GetModelID()<<"\t"<<result.GetTrainingT(0)<<"\n";
}
float TFinalReport::GetBestMOR(string& id) {
	float bestValue = 0.0;
	
	for (int i=EvaluationMethod::minModelSize; i<(int)reports.size(); i++) {
		string curID;

		float curValue = GetBestMOR(i, curID);
		if (curValue > bestValue) {
			id=curID;
			bestValue=curValue;
		}
	}
	return bestValue;
}
float TFinalReport::GetBestMOR(int order, string& id) {
	FoldStatistic bestModel = reports[order].GetBestModel(0);
	id = bestModel.model.GetModelID();
	
	//If we want the biased one, we can just return this
	//return bestModel.AvgMOR();
	return reports[order].GetAvgMOR(0);
}

void TFinalReport::GenerateVerboseReport(std::ostream& os, int depthOfReport, EvaluationMethod *eval) {
	int orderSize = reports.size();
	for (int i=0; i<depthOfReport; i++) 
		for (int o=0; o<orderSize;o++) 
			reports[0].GenerateDetailedReport(i, os, eval);
}

void TFinalReport::GenerateSummaryReport(std::ostream& os, Distribution::PTestDistribution* dist) {
	for (uint i=0; i<reports.size(); i++) 
		reports[i].GenerateSummaryReport(os, dist);
}	

void TFinalReport::GenerateDetailedReport(int idx, std::ostream& os, EvaluationMethod* eval) {
	for (uint i=0; i<reports.size(); i++) 
		reports[i].GenerateDetailedReport(idx, os, eval);
}

}

}
