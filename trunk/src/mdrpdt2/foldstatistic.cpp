//
// C++ Implementation: foldstatistic
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "foldstatistic.h"

namespace MdrPDT {

namespace Evaluation {

float FoldEvaluation::reportThreshold = 0.0;

void FoldEvaluation::AddFold(PdtModel& report, int foldID) {
	stats[report.GetModelID()].AddFold(report, foldID);
	modelIDs.push_back(report.GetModelID());
	
	if ((int)modelIDs.size()-1 != foldID)
		cout<<"AddFold() "<<modelIDs.size()<<" != "<<foldID<<"\n";
	assert((int)modelIDs.size()-1 == foldID);			///<Just making sure we don't re-add models
}

void FoldEvaluation::ReportHeader(std::ostream& os, int foldCount, Distribution::PTestDistribution *dist) {
	int width = 10;

	int labelWidth = 2 * PdtModel::MaxLabelLength + 4;	

	os<<"\n\n";
	///////////////////////// 1rst Line
	os<<setw(8)<<" ";
	if (foldCount > 1)
		os<<setw(width)<<" ";
	os<<setw(labelWidth)<<" ";
	os<<setw(width)<<"MDR-PDT";
	if (foldCount>1)
		os<<setw(width)<<"MDR-PDT";
	os<<setw(width)<<"Matched";
	os<<"\n";
	///////////////////////// 2nd Line
	os<<setw(8)<<"Model";
	if (foldCount > 1)
		os<<setw(width)<<"XV";
	os<<setw(labelWidth)<<" ";
	if (foldCount>1)
		os<<setw(width)<<"Train."<<setw(width)<<"Test";
	else
		os<<setw(width)<<"T";
	os<<setw(width)<<"Odds";
	if (dist)
		os<<setw(width)<<"MOR";
	os<<"\n";
	///////////////////////// 3rd line		
	os<<setw(8)<<"Rank";
	if (foldCount > 1)
		os<<setw(width)<<"Cons.";

	os<<setw(labelWidth)<<"Model";
	os<<setw(width)<<"Stat.";
	if (foldCount > 1)
		os<<setw(width)<<"Stat.";

	os<<setw(width)<<"Ratio";
	if (dist)
		os<<setw(width)<<"P-Value";
	os<<"\n-------------------------------------------------------------------------------\n";
}

FoldStatistic FoldEvaluation::GetBestModel() {
	map<string, FoldStatistic>::iterator itr = stats.begin();
	map<string, FoldStatistic>::iterator end = stats.end();
	
	FoldStatistic bestModel = itr->second;
	while (++itr != end) {
		if (bestModel < itr->second)
			bestModel = itr->second;
	}
	return bestModel;
}

PdtModel FoldEvaluation::ReportBestModel(std::ostream& os, int rank, int foldCount, Distribution::PTestDistribution *distribution) {
	int ptestSigDigits = 3;
	if (distribution)
		ptestSigDigits = distribution->GetSignificantDigits();
	int width = 10;
	int labelWidth = 2 * PdtModel::MaxLabelLength + 4;	
	FoldStatistic bestModel = GetBestModel();
	float avgMOR = bestModel.AvgMOR();
	float sumMOR = 0.0;
	stringstream ss;
	if (foldCount > 1) {
		for (int f=0; f<foldCount; f++) {
			PdtModel &model = stats[modelIDs[f]].model;
			ss<<setw(8)<<" "<<setw(width-1)<<" "<<" ";		//Rank is ignored for the individual slices / CVC
			ss<<setw(labelWidth)<<model.GetPrintableLabel();
			ss<<setw(width)<<setprecision(3)<<model.GetTrainingT(f);
			ss<<setw(width)<<setprecision(3)<<model.GetTestingT(f);
			MatchedOddsRatio mor = model.EvaluateMOR(f);
			ss<<setw(width)<<setprecision(3)<<mor.GetRatio();
			sumMOR += mor.GetRatio();
//The biased stuff seems to have a memory error associated with it. I need to fix that...but later
//			os<<setw(width)<<setprecision(3)<<model.EvaluateBiasedMOR(f).GetRatio();
			ss<<"\n";
		}
		ss<<"-------------------------------------------------------------------------------\n";
	} else {
		sumMOR = bestModel.model.EvaluateMOR().GetRatio();
	}
	avgMOR = sumMOR / (float)foldCount;
	float pValue = distribution->Evaluate(bestModel.model.ModelSize()-1,avgMOR);
	if (reportThreshold > 0.0 && pValue > reportThreshold)
		return bestModel.model;
	os<<ss.str();
	os<<setw(8)<<rank;
	if (foldCount > 1)
		os<<setw(width-1)<<bestModel.GetCount()<<" ";
	os<<setw(labelWidth)<<bestModel.model.GetPrintableLabel();
	os<<setw(width)<<setprecision(3)<<bestModel.AvgTrainingT();
	if (foldCount >1) 
		os<<setw(width)<<setprecision(3)<<bestModel.AvgTestingT();
	os<<setw(width)<<setprecision(3)<<sumMOR / (float)foldCount;
	//os<<setw(width)<<setprecision(3)<<bestModel.AvgMOR();
	if (distribution) 
		os<<setw(width)<<setprecision(ptestSigDigits)<<pValue;
	os<<"\n\n";
	return bestModel.model;	
}

FoldStatistic::FoldStatistic() : foldCount(0), cvc(0.0), avgMOR(0.0), sumTrainingT(0.0), sumTestingT(0.0) { }

FoldStatistic::FoldStatistic(int xvCount) : foldCount(xvCount), cvc(0.0), avgMOR(0.0), sumTrainingT(0.0), sumTestingT(0.0) { }

int FoldStatistic::GetXVConsitency() { return (int)cvc; }

void FoldStatistic::AddFold(PdtModel& model, int foldID) {
	this->model=model;
	foldIDs.push_back(foldID);
	float training, testing;
	model.EvaluateFold(training, testing, foldID);
	sumTrainingT+=training;
	sumTestingT+=testing;
	cvc+=1.0;
}


float FoldStatistic::AvgTrainingT() const {	return sumTrainingT / cvc;	}


float FoldStatistic::AvgTestingT() const { return sumTestingT / cvc; }

float FoldStatistic::AvgBiasedMOR() {
//	if (avgMOR == 0.0) {
		vector<int>::iterator itr = foldIDs.begin();
		vector<int>::iterator end = foldIDs.end();
		float avgMOR = 0.0, sumMOR = 0.0;
		while (itr!=end) {
			sumMOR+=model.EvaluateMOR(*itr).GetRatio();
			itr++;
		}
		avgMOR = sumMOR/((float)(foldIDs.size()));
//	}
	return avgMOR;
}

float FoldStatistic::SumMOR() {
	int foldCount = foldIDs.size();
	float sumMOR = 0.0;
	for (int i=0; i<foldCount; i++) {
		MatchedOddsRatio mor;
		if (foldCount>1)
			sumMOR+=model.EvaluateMOR(foldIDs[i]).GetRatio();
		else
			sumMOR+=model.EvaluateMOR(foldIDs[0]).GetRatio();
	}
	return sumMOR;
}
//This requires a bit more work, so we'll only do it for the cases we care about
float FoldStatistic::AvgMOR() {
	int foldCount = foldIDs.size();
	float avgMOR = 0.0, sumMOR = 0.0;
	if (this->avgMOR == 0.0) {
		for (int i=0; i<foldCount; i++) {
			MatchedOddsRatio mor;
			if (foldCount>1)
				sumMOR+=model.EvaluateMOR(foldIDs[i]).GetRatio();
			else
				sumMOR+=model.EvaluateMOR(foldIDs[0]).GetRatio();
		}
		avgMOR = sumMOR/((float)(foldCount));
	}
	return avgMOR;
}
	
void FoldStatistic::SetAvgMOR(float avgMOR) {
	this->avgMOR = avgMOR;
}

bool FoldStatistic::operator<(FoldStatistic& other) {
	if (cvc == other.cvc)
		return AvgTestingT() < other.AvgTestingT();
	return cvc < other.cvc;

	//return count < other.count || AvgTestingT() < other.AvgTestingT() || AvgMOR() < other.AvgMOR();
}


int FoldStatistic::GetCount() { return (int)cvc; }



}

}

