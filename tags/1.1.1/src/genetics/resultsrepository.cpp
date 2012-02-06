//
// C++ Implementation: resultsrepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "resultsrepository.h"
#include <iomanip>

namespace Genetics {

namespace Reporting {

ResultsRepository::~ResultsRepository()
{
	delete[] results;
}




ModelStatistics ResultsRepository::GetModelStatistics(uint order, uint position) {
	ModelStatistics stats;

	assert(order < maxModelSize);
	for (uint i=0; i<foldCount; i++) {
		uint idx = order * foldCount + i;
		assert(position < results[idx].size());
		stats.Append(i, results[idx][position]);		
	}
	return stats;
}





float ResultsRepository::GetAverageTesting(uint modelSize, uint position, float &oddsRatio, string& snpLabel) {
	float sumT = 0.0;
	map<string, TestHolder> frequencies;
	for (uint i=0; i<foldCount; i++) {
		uint idx = modelSize * foldCount + i;

		FoldStatistics &stat = results[idx][position];
		frequencies[stat.label].AddTest(stat);
		sumT += stat.testingFitness;
	}

	vector<TestHolder> winners;					///<We'll stash each of the winners to find the real winner
	map<string, TestHolder>::iterator itr = frequencies.begin();
	map<string, TestHolder>::iterator end = frequencies.end();
	
	for (; itr!=end; itr++) 
		winners.push_back(itr->second);
		//winners[i++] = itr->second;

	sort(winners.begin(), winners.end());
	//float value = winners[winners.size() - 1].AverageFitness();
	float value = sumT / (float)foldCount;
	snpLabel=winners[winners.size() - 1].snpLabel;	
	oddsRatio = winners[winners.size() - 1].oddsRatio;
	return value;

}


/*
	 	   		   Avg      Avg          Matched   
Model	Fold	    T      Class.          Odds     p
Size	  # 	  Stat.    Error          Ratio   Value
----------------------------------------------------------
 

*Fold # is optional based on foldCount > 1
*/
void ResultsRepository::WriteHeader(ostream *os) {
	uint width = 12;
	
	uint firstJump = 20;
	
	if (foldCount > 1)
		firstJump+=24;
	*os<<setw(firstJump + width)<<"MDR-PDT";
	
	if (foldCount > 1)
		*os<<setw(width)<<""<<setw(width)<<"MDR-PDT";

	*os<<setw(width)<<""<<setw(width)<<"Matched"<<endl;

	*os<<setw(8)<<"Model";
	
	if (foldCount > 1)
		*os<<setw(8)<<"Model"<<setw(8)<<"XV"<<setw(8)<<"Fold";
	*os<<setw(width)<<""<<setw(width)<<"Train."<<setw(width)<<"Class.";
	if (foldCount > 1)
		*os<<setw(width)<<"Test"<<setw(width)<<"Pred.";
	*os<<setw(width)<<"Odds ";
	if (fitnessDist)
		*os<<setw(width+6)<<fitnessDist->GetLabel();
	if (orDist)
		*os<<setw(width+6)<<orDist->GetLabel();
	*os<<"\n";

	*os<<setw(8)<<"Rank";
	if (foldCount > 1)
		*os<<setw(8)<<"Size"<<setw(8)<<"Cons."<<setw(8)<<"# ";

	*os<<setw(width)<<"Model";
	if (foldCount > 1)
		*os<<setw(width)<<"Stat."<<setw(width)<<"Error";

	*os<<setw(width)<<"Stat."<<setw(width)<<"Error"<<setw(width)<<"Ratio";
	if (fitnessDist)
		*os<<setw(width + 5)<<"p Value";
	if (orDist)
		*os<<setw(width + 5)<<"p Value";

	*os<<endl;

	*os<<"------------------------------------------------------------";
	if (foldCount > 1)
		*os<<"-----------------------------------------";
	if (fitnessDist)
		*os<<"---------------------";
	if (orDist)
		*os<<"---------------------";
	*os<<"\n";
}

void ResultsRepository::BakeResultModels(uint start, uint stop, uint reportSize, SnpRepository *repo, SnpEvaluation *eval) {

//	summaryReport.clear();
/*	for (uint entry = 0; entry < reportSize; entry++) {
		for (uint p=start; p<stop; p++){
			ModelStatistics stats = GetModelStatistics(p, entry);
			stats.GatherTestStatistics(repo, eval);
//			summaryReport.push_back(stats);
		}
	}*/
}


void ResultsRepository::GenerateReport(ostream *os, uint start, uint stop, uint reportSize, SnpRepository *repo, SnpEvaluation *eval) {
	uint width = 12;

	if (reportSize == 0) 
		reportSize = 5;

	uint idx = start * foldCount;

	if (reportSize > results[idx].size())
		reportSize = results[idx].size();

	if (reportSize == 0)
		*os<<"No data to summarize\n";
	
	WriteHeader(os);
	Sort();


	for (uint entry = 0; entry < reportSize; entry++) {
		for (uint p=start; p<stop; p++){
			ModelStatistics stats = GetModelStatistics(p, entry);

			//If the fold count is 1, we don't really need to duplicate the information
			if (foldCount > 1)
				stats.GatherTestStatistics(repo, eval, os);
			else
				stats.GatherTestStatistics(repo, eval, NULL);
			//stats.GatherTestStatistics(repo, eval);
			float testStat = stats.GetAvgTesting();
			*os<<setw(8)<<entry+1;

			if (foldCount > 1)
				*os<<setw(16)<<stats.xvConsistency<<setw(8)<<"";
			*os<<setw(width)<<stats.GetLabel();
			if (foldCount > 1) {
				*os<<setw(width)<<setprecision(3)<<stats.GetAvgTraining();
				*os<<setw(width)<<setprecision(3)<<stats.GetAvgClassificationError();
			}
			*os<<setw(width)<<setprecision(3)<<testStat;
			*os<<setw(width)<<setprecision(3)<<stats.GetAvgPredictionError();
			*os<<setw(width)<<setprecision(3)<<stats.GetOddsRatio();
			if (fitnessDist)
				*os<<setw(width)<<"     p < "<<(fitnessDist->GetPValue(testStat, p));
			if (orDist)
				*os<<setw(width)<<"     p < "<<(orDist->GetPValue(stats.GetOddsRatio(), p));
			if (peDist)
				*os<<setw(width)<<"     p < "<<(peDist->GetPValue(stats.GetAvgPredictionError(), p));
			*os<<"\n";
		}
	}
	
}

void ResultsRepository::GenerateReport(ostream* os, uint reportSize, uint start, uint stop) {
	assert(stop < maxModelSize);
	Sort();

	uint width = 12;
	WriteHeader(os);

	if (reportSize == 0) 
		reportSize = 5;
	uint idx = start * foldCount;

	if (reportSize > results[idx].size())
		reportSize = results[idx].size();

	if (reportSize == 0)
		*os<<"No data to summarize\n";

	for (uint entry = 0; entry < reportSize; entry++) {
		for (uint p=start; p<stop; p++){
			string label;
			float oddsRatio = 0.0;
			float avgClsError = 0.0;
			float avgPredError = 0.0;
	
	
			for (uint i=0; i<foldCount; i++) {
				uint idx = p * maxModelSize + i;
				cout<<"I need to print out the break down for each fold\n";
				avgClsError+=results[idx][entry].clsError;
				avgPredError+=results[idx][entry].predError;
			}
			float avgFitness = GetAverageTesting(p, entry, oddsRatio, label);

			*os<<p+1<<"\t";
			*os<<setw(width)<<label;
			*os<<setw(width)<<avgFitness;
			*os<<setw(width)<<setprecision(3)<<avgClsError/(float)foldCount;
			*os<<setw(width)<<setprecision(3)<<avgPredError/(float)foldCount;
			*os<<setw(width)<<setprecision(3)<<oddsRatio;
			if (fitnessDist)
				*os<<"     p < "<<(fitnessDist->GetPValue(avgFitness, p))<<"\n";
			if (orDist)
				*os<<"     p < "<<(orDist->GetPValue(oddsRatio, p))<<"\n";
		}
	}

}










}

}
