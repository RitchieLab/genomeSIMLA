//
// C++ Implementation: evalbalancedaccuracy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "evalbalancedaccuracypdt.h"
#include "matchedoddsratio.h"
namespace MDR {

bool EvalBalancedAccuracyPDT::Verbose = false;

void EvalBalancedAccuracyPDT::SummarizeResults(ostream *os) {
	uint width = 12;
	*os<<"Model\n";
	*os<<"Size\t";

	*os<<"Fold\t";

	*os<<setw(width)<<"Model"<<setw(width)<<"T Statstic";
	if (fitnessDist)
		*os<<setw(width)<<"p Value";
	*os<<endl;

	for (uint p=comboStart-1; p<comboEnd; p++){
		for (uint i=0; i<foldCount; i++) {
			if (fitnessDist) {
				*os<<p+1<<"\t";

				*os<<i+1<<"\t";

				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].topStatistics[p].trainingFitness;
				if (foldCount > 0)
					*os<<setw(width)<<results[i].topStatistics[p].testingFitness;
				*os<<setw(width)<<results[i].topStatistics[p].clsError;
				*os<<setw(width)<<results[i].topStatistics[p].predError;
				*os<<setw(width)<<results[i].topStatistics[p].oddsRatio;
				*os<<"     p < "<<fitnessDist->GetPValue(results[i].topStatistics[p].trainingFitness, results[i].topModel[p]->GetLabelCount())<<"\n";
			}
			else {
				*os<<p+1<<"\t";

				*os<<i+1<<"\t";

				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].topStatistics[p].trainingFitness;
				if (foldCount > 0)
					*os<<setw(width)<<results[i].topStatistics[p].testingFitness;
				*os<<setw(width)<<results[i].topStatistics[p].clsError;
				*os<<setw(width)<<results[i].topStatistics[p].predError;
				*os<<setw(width)<<results[i].topStatistics[p].oddsRatio;
			}
		}
	}
}


void EvalBalancedAccuracyPDT::BakeResultModels(uint reportSize, SnpRepository *repo) {
	if (overallResults)
		overallResults->BakeResultModels(comboStart-1, comboEnd, reportSize, repo, this);

}

void EvalBalancedAccuracyPDT::ReportResults(ostream *os, uint reportSize, SnpRepository *repo) {
	if (overallResults) {
		overallResults->GenerateReport(os, comboStart-1, comboEnd, reportSize, repo, this);
	}
}

void EvalBalancedAccuracyPDT::ReportResults(ostream *os) {
	*os<<"Top Model of the PDT Search:\n";
	uint width = 12;

	*os<<"    \t"<<setw(width)<<""<<setw(width)<<"Avg "<<setw(width)<<"Avg" <<setw(width)<<"Avg "<<setw(width)<<"Matched"<<endl;
	*os<<"Model\t"<<setw(width)<<" "<<setw(width)<<"MDR-PDT"<<setw(width)<<"Class."<<setw(width)<<"Pred."<<setw(width)<<"Odds ";
	if (fitnessDist)
		*os<<setw(width)<<"Avg"<<"\n";
	*os<<"Size\t";
	*os<<setw(width)<<"Model"<<setw(width)<<"Stat."<<setw(width)<<"Error"<<setw(width)<<"Error"<<setw(width)<<"Ratio";
	if (fitnessDist)
		*os<<setw(width)<<"p Value";
	*os<<endl;
	for (uint p=comboStart-1; p<comboEnd; p++){
		SnpAligned *snp 	= NULL;
		float oddsRatio 	= 0.0;
//		float fitness = 0.0;
		float avgFitness = GetAverageTesting(p, oddsRatio, snp);
		float avgClsError 	= 0.0;
		float avgPredError 	= 0.0;
		float avgMO			= 0.0;
		

		*os<<p+1<<"\t";
		*os<<setw(width)<<snp->GetLabel();
		*os<<setw(width)<<avgFitness;

		for (uint i=0; i<foldCount; i++) {
			avgClsError+=results[i].topStatistics[p].clsError;
			avgPredError+=results[i].topStatistics[p].predError;
			avgMO+=results[i].topStatistics[p].oddsRatio;
		}
		*os<<setw(width)<<setprecision(3)<<avgClsError/(float)foldCount;
		*os<<setw(width)<<setprecision(3)<<avgPredError/(float)foldCount;
		*os<<setw(width)<<setprecision(3)<<avgMO/(float)foldCount;
		if (fitnessDist)
 			*os<<"     p < "<<fitnessDist->GetPValue(avgFitness, p)<<"\n";
	}
}

void EvalBalancedAccuracyPDT::DisplayAnalyses(bool writeDetails) {
	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();

	stringstream details;

	*os<<"\n"<<setw(55)<<right<<" Matched\n";
	if (foldCount > 1)
		*os<<setw(5)<<"Fold";

	*os<<setw(20)<<right<<"Model ";
	*os<<setw(10)<<right<<"Class.";
	if (foldCount > 1)  {
		*os<<setw(20)<<right<<"MDR-PDT Stat ";
		*os<<setw(10)<<right<<"Pred. ";
	}
	else
		*os<<setw(10)<<right<<"MDR-PDT ";

	*os<<setw(10)<<right<<"Odds";

	*os<<endl;

	if (foldCount > 1)
		*os<<setw(5)<<"  #";

	*os<<setw(20)<<right<<"  ID ";
	*os<<setw(10)<<right<<"Error";
	if (foldCount > 1) {
		*os<<setw(10)<<right<<"Training";
		*os<<setw(10)<<right<<"Testing";
		*os<<setw(10)<<right<<"Error";
	}
	else {
		*os<<setw(10)<<right<<"Stat";
	}

	*os<<setw(10)<<right<<"Ratio";

	*os<<endl;
	*os<<"--------------------------------------------------------------------\n";
	
	

	//Iterate over each of the status objects insde the vector
	for (uint i=0; i<foldCount; i++) {

		details<<"\nModel Details - Training Set #"<<i<<" of "<<foldCount<<"\n";
		details<<"------------------------------------------------------\n";

		for (uint p=comboStart-1; p<comboEnd; p++){
			Genetics::Evaluation::FoldStatistics &stats = results[i].topStatistics[p];
			stats = EvaluateVerbose(i, results[i].topModel[p], details);
			//results[i].topStatistics[p] = newStats;

		}
		
	}
	
	if (writeDetails && !isPTest) {
		details<<"\t\t*Indicates models that are determined to be High Risk\n";
		details<<"\n\t\t Heterozygote alleles are not necessarily ordered the\n";
		details<<"\t\t same as they were found in the original data\n";
		*os<<details.str()<<"\n";

	}
}

//Eventually, this will probably just set up the scores
void EvalBalancedAccuracyPDT::PostEvaluation() {
	//Iterate over each of the status objects insde the vector
/*	for (uint i=0; i<foldCount; i++) {
		for (uint p=comboStart-1; p<comboEnd; p++){
			if (isPTest) {
				//cout<<results[i].highTStatistic.size();
				if (results[i].topStatistics[p].trainingFitness > 0.0001)
					distribution->AppendTest(results[i].topStatistics[p].trainingFitness, results[i].topModel[p]);	
			}
		}
		
	}*/
}


FoldStatistics EvalBalancedAccuracyPDT::EvaluateVerbose(uint fold, SnpAligned *snp, stringstream& details) {
	uint totalAffected = 0;
	uint totalUnaffected = 0;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	//uint cellCount=0;
	Genetics::Evaluation::FoldStatistics bestStatistics;
	int  sumD = 0;
	int sumD2 =0;
	int localD = 0;						//This is the D that represents the local slice only
	int localD2 = 0;						//This represents the D2 for the local slice
	
	if (snp == NULL)
		return bestStatistics;
	
	uint gtCount = snp->CountGenotypes();
	//uint highRiskCells[gtCount];
	BitSetType highRisk;
	highRisk.resize(gtCount);

	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();
	string tag;
	sumD2 = 0;
	sumD  = 0;
	//cellCount = 0;

	PdtFold &validation = folds[fold];
	CaseControlStatus training = validation.GetOverallStatus();
	CaseControlStatus testing  = validation.GetOverallStatus();

	if (foldCount > 1)
		training = validation.GetTrainingSet();

		
	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);

	//totalAffected = (curGenotype & training.affected).count();
	//totalUnaffected = (curGenotype & training.unaffected).count();

	float totalIndividuals = (float)training.affected.size();

	//cout<<"Training Mask: "<<training.total<<"\n";
	details<<setw(20)<<right<<"\nGenotype"<<setw(18)<<right<<"Individuals      ";
	stringstream gtLabel;
	for (uint n=0; n<snp->GetLabelCount(); n++) 
		gtLabel<<setw(4)<<right<<snp->GetLabel(n)<<" ";

	details<<"\n"<<setw(20)<<right<<gtLabel.str();
	details<<"  "<<setw(6)<<"A"<<setw(6)<<"U";
	details<<setw(6)<<"TOT"<<setw(8)<<"Ratio"<<endl;      
	details<<"Total Individuals in Training Set: "<<training.total.count()<<"\n";
	//details<<"Total Individuals in Testing Set:  "<<testing.affected.count()<<"\n";      
	details<<"    ------------------------------------------------\n";


	highRisk.reset();

	for (uint gtID=1; gtID<gtCount; gtID++) {

		curGenotype=snp->GetGenotype(gtID);
		//Determine if the cells are high risk
		localAffected = (curGenotype & training.affected).count();
		localUnaffected = (curGenotype & training.unaffected).count();

		if (localUnaffected >0) {
			localRatio = (float)localAffected/(float)localUnaffected;		
		}
		else
			localRatio=localAffected;
		totalAffected += localAffected;
		totalUnaffected += localUnaffected;

		if (localRatio >= 1.0) {
			tag=" *";
			sumD+=(localAffected - localUnaffected);
			localD+=((curGenotype & testing.affected).count() - (curGenotype & testing.unaffected).count());
			highRisk[gtID] = true;
		}		
		else {
			tag="  ";
		}
		details<<setw(20)<<right<<snp->GetGenotypeLabel(gtID)<<tag<<setw(6)<<(uint)localAffected;
		details<<setw(6)<<(uint)localUnaffected<<setw(6)<<((uint)localAffected+(uint)localUnaffected);
		details<<setw(8)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(2)<<localRatio<<endl;
	}
	snp->SetHrCells( highRisk );
	details<<endl;
	details<<"Model Details [ "<<snp->GetLabel()<<" ]\n";
	details<<"\t    Affected: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<((float)totalAffected/totalIndividuals)*100.0 <<"%\t";
	details<<"\tUnaffected: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<((float)totalUnaffected/totalIndividuals)*100.0<<"%\n";
	details<<"\tMissing Data: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<(100.0 - ((float)(totalAffected+totalUnaffected)/(float)curGenotype.size())*100.0)<<"% -- "<<((float)(~(snp->GetGenotype(0))).count()/(float)curGenotype.size() * 100.0)<<"\n";

	//Risk Cells are built, let's build the score for these pieces
	for (uint j=0; j<foldCount; j++) {
		//Let's skip the one that is the validation set
		if (fold != j || foldCount == 1) {
			sumD2+=folds[j].CalculateD2(snp);
		}
	}

	localD2 = folds[fold].CalculateD2(snp);

	validation.results.CalcTrainingT(sumD, sumD2);
	validation.CalculateAccuracy(snp);
	validation.results.CalcTestingT(localD, localD2);
	validation.ResetCache();
	

	if (foldCount > 1)
		*os<<setw(5)<<fold;

	float classificationError = validation.results.ClassificationError();
	float predictionError = validation.results.PredictionError();

	if (foldCount == 1)
		classificationError = predictionError;

	float matchedOddsRatio = 0.0;


	*os<<setw(19)<<right<<snp->GetLabel()<<" ";
	*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<right<<classificationError<<"%";
	*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<right<<validation.results.trainingT;
	if (foldCount > 1) {
		*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<right<<validation.results.testingT;
		*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<right<<predictionError<<"%";
	}


	if (sibshipStatus.size() > 0)
		matchedOddsRatio = FindMatchedOddsRatio(snp, sibshipStatus[0], training);
	*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3)<<right<<matchedOddsRatio;



	*os<<endl;

	bestStatistics.testingFitness = validation.results.testingT;
	bestStatistics.trainingFitness = validation.results.trainingT;
	bestStatistics.clsError = classificationError;
	bestStatistics.predError = predictionError;
	bestStatistics.oddsRatio = matchedOddsRatio;
	

	//if (validation.results.trainingT > bestT)
	//	bestT = validation.results.trainingT;

	return bestStatistics;
}

FoldStatistics EvalBalancedAccuracyPDT::GetFoldStatistics(uint fold, SnpAligned *snp) {
	uint totalAffected;
	uint totalUnaffected;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	Genetics::Evaluation::FoldStatistics bestStatistics;
	int sumD = 0;
	int sumD2 =0;
	int localD = 0;						//This is the D that represents the local slice only
	int localD2 = 0;						//This represents the D2 for the local slice
	
	if (snp == NULL)
		return bestStatistics;
	uint gtCount = snp->CountGenotypes();
	BitSetType highRisk(gtCount, 0);

	string tag;
	sumD2 = 0;
	sumD  = 0;

	PdtFold &validation = folds[fold];
	CaseControlStatus training = validation.GetOverallStatus();
	CaseControlStatus testing  = validation.GetOverallStatus();
	
	if (foldCount > 1)
		training = validation.GetTrainingSet();
		
	//Figure out what the overall ratio is
	BitSetType curGenotype = ~snp->GetGenotype(0);

	totalAffected = (curGenotype & training.affected).count();
	totalUnaffected = (curGenotype & training.unaffected).count();

//	float totalIndividuals = (float)training.affected.size();
	highRisk.reset();

	for (uint gtID=1; gtID<gtCount; gtID++) {
		curGenotype=snp->GetGenotype(gtID);
		//Determine if the cells are high risk
		localAffected = (curGenotype & training.affected).count();
		localUnaffected = (curGenotype & training.unaffected).count();

		if (localUnaffected >0) {
			localRatio = (float)localAffected/(float)localUnaffected;		
		}
		else
			localRatio=localAffected;
		
		if (localRatio >= 1.0) {
			sumD+=(localAffected - localUnaffected);
			localD+=((curGenotype & testing.affected).count() - (curGenotype & testing.unaffected).count());
			highRisk[gtID] = true;
		}		
	}

	snp->SetHrCells(highRisk);

	//Risk Cells are built, let's build the score for these pieces
	for (uint j=0; j<foldCount; j++) {
		//Let's skip the one that is the validation set
		if (fold != j || foldCount == 1) {
			sumD2+=folds[j].CalculateD2(snp);
		}
	}
	localD2 = folds[fold].CalculateD2(snp);
	validation.results.CalcTrainingT(sumD, sumD2);
	validation.CalculateAccuracy(snp);
	validation.results.CalcTestingT(localD, localD2);
	validation.ResetCache();
	float classificationError = validation.results.ClassificationError();
	float predictionError = validation.results.PredictionError();
	float matchedOddsRatio = 0.0;



	if (sibshipStatus.size() > 0) 
		matchedOddsRatio = FindMatchedOddsRatio(snp, sibshipStatus[0], training);


	bestStatistics.testingFitness = validation.results.testingT;
	bestStatistics.trainingFitness = validation.results.trainingT;
	bestStatistics.clsError = classificationError;
	bestStatistics.predError = predictionError;
	bestStatistics.oddsRatio = matchedOddsRatio;
	bestStatistics.label = snp->GetLabel();
	bestStatistics.order = snp->GetLabelCount();
	//if (validation.results.trainingT > bestT)
	//	bestT = validation.results.trainingT;

	return bestStatistics;
}
/*
FoldStatistics EvalBalancedAccuracyPDT::EvaluateFold(uint fold, SnpAligned *snp) {
	uint totalAffected;
	uint totalUnaffected;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	Genetics::Evaluation::FoldStatistics bestStatistics;
	int  sumD = 0;
	int sumD2 =0;
	int localD = 0;						//This is the D that represents the local slice only
	int localD2 = 0;						//This represents the D2 for the local slice
	
	uint gtCount = snp->CountGenotypes();
	BitSetType highRisk(gtCount);

	string tag;
	sumD2 = 0;
	sumD  = 0;
	cellCount = 0;

	PdtFold &validation = folds[fold];
	CaseControlStatus training = validation.GetOverallStatus();
	//CaseControlStatus testing  = validation.GetOverallStatus();
	
	if (foldCount > 1)
		training = validation.GetTrainingSet();
		
	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);

	totalAffected = (curGenotype & training.affected).count();
	totalUnaffected = (curGenotype & training.unaffected).count();

//	float totalIndividuals = (float)training.affected.size();
	
	for (uint gtID=1; gtID<gtCount; gtID++) {

		curGenotype=snp->GetGenotype(gtID);
		//Determine if the cells are high risk
		localAffected = (curGenotype & training.affected).count();
		localUnaffected = (curGenotype & training.unaffected).count();

		if (localUnaffected >0) {
			localRatio = (float)localAffected/(float)localUnaffected;		
		}
		else
			localRatio=localAffected;
		
		if (localRatio >= 1.0) {
			//sumD+=(localAffected - localUnaffected);
			//localD+=((curGenotype & testing.affected).count() - (curGenotype & testing.unaffected).count());
			highRisk[gtId] = true;
		}		
	}

	//Risk Cells are built, let's build the score for these pieces
	for (uint j=0; j<foldCount; j++) {
		//Let's skip the one that is the validation set
		if (fold != j || foldCount == 1) {
			sumD2+=folds[j].CalculateD2(snp);
		}
	}

#ifdef USE_DOPTIMIZATION
	if (!isPTest || qdMax[snp->GetLabelCount() - 1].Evaluate(sumD)) 
#endif
	{

		localD2 = folds[fold].CalculateD2(snp);
	
		validation.results.CalcTrainingT(sumD, sumD2);
		validation.CalculateAccuracy(snp);
		validation.results.CalcTestingT(localD, localD2);
		validation.ResetCache();
		
		float classificationError = validation.results.ClassificationError();
		float predictionError = validation.results.PredictionError();
		float matchedOddsRatio = 0.0;



#ifdef CROSS_VALIDATION
		if (sibshipStatus[0] > 0)
			matchedOddsRatio = FindMatchedOddsRatio(snp, sibshipStatus[0]);


#endif

		bestStatistics.testingFitness = validation.results.testingT;
		bestStatistics.trainingFitness = validation.results.trainingT;
		bestStatistics.clsError = classificationError;
		bestStatistics.predError = predictionError;
		bestStatistics.oddsRatio = matchedOddsRatio;
	

	//if (validation.results.trainingT > bestT)
	//	bestT = validation.results.trainingT;
	}
	return bestStatistics;
}*/

BitSetType EvalBalancedAccuracyPDT::GetHighRiskCells( SnpAligned * snp) {
	CaseControlStatus ind = folds[0].GetOverallStatus();

	for (uint i=1; i<foldCount; i++) {
		CaseControlStatus fold = folds[i].GetOverallStatus();
		ind.affected 	|= fold.affected;
		ind.unaffected 	|= fold.unaffected;
		ind.total  		|= fold.total;
	}

	return GetHighRiskCells(snp, ind);
}
	
/**
 * @brief Quickly walks through and determines the high risk cells for a given set of individuals
 */
BitSetType EvalBalancedAccuracyPDT::GetHighRiskCells(SnpAligned *snp, CaseControlStatus &individuals) {
	uint localAffected;
	uint localUnaffected;
	float localRatio;
	
	//cout<<snp->GetLabel()<<":Affected   = "<<individuals.affected<<"\n";
	//cout<<snp->GetLabel()<<":Unaffected = "<<individuals.unaffected<<"\n";


	uint gtCount = snp->CountGenotypes();
	BitSetType highRisk(gtCount);

	//Figure out what the overall ratio is
	BitSetType curGenotype;

	//CaseControlStatus trainingSet = 	//individuals
	snp->ResetEvaluations();

	//curGenotype = snp->GetGenotype(0);
	highRisk.reset();
	for (uint gtID=1; gtID<gtCount; gtID++) {
		curGenotype=snp->GetGenotype(gtID);
		//Determine if the cells are high risk
		localAffected = (curGenotype & individuals.affected).count();
		localUnaffected = (curGenotype & individuals.unaffected).count();

		if (localUnaffected >0) {
			localRatio = (float)localAffected/(float)localUnaffected;		
		}
		else
			localRatio=localAffected;
		
		if (localRatio >= 1.0) {
			highRisk[gtID] = true;
		}
				
	}

	//This is debatable as to whether we should do this
	snp->SetHrCells(highRisk);
	
	return highRisk;
}

bool EvalBalancedAccuracyPDT::EvaluateModel(SnpAligned *snp, uint position) {
	uint  sumD = 0;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	uint sumD2 =0;
	
	uint gtCount = snp->CountGenotypes();
	BitSetType highRisk(gtCount);
	float bestT = 0.0;

//	static int lastmodel = -1;

	//Figure out what the overall ratio is
	BitSetType curGenotype;

	CaseControlStatus trainingSet = folds[0].GetOverallStatus();
	snp->ResetEvaluations();
	for (uint i=0; i<foldCount; i++) {
		sumD2 = 0;
		sumD  = 0;

		PdtFold &validation = folds[i];
		if (foldCount > 1) 
			trainingSet = validation.GetTrainingSet();

/*		if (position == 1) {
			cout<<"EvaluateModel("<<snp->GetLabel()<<") : ";
			cout<<"    Affected: "<<trainingSet.affected<<"\n";
			lastmodel = snp->GetLabelCount();
		}
*/
		//curGenotype = snp->GetGenotype(0);
		highRisk.reset();
		for (uint gtID=1; gtID<gtCount; gtID++) {
			curGenotype=snp->GetGenotype(gtID);
			//Determine if the cells are high risk
			localAffected = (curGenotype & trainingSet.affected).count();
			localUnaffected = (curGenotype & trainingSet.unaffected).count();

			if (localUnaffected >0) {
				localRatio = (float)localAffected/(float)localUnaffected;		
			}
			else
				localRatio=localAffected;
			
			if (localRatio >= 1.0) {
				sumD+=(localAffected - localUnaffected);
				highRisk[gtID] = true;
			}
				
		}

		snp->SetHrCells(highRisk);
		//cout<<"-->"<<snp->GetLabel()<<" : "<<highRisk<<"\n";
#ifdef USE_DOPTIMIZATION
		if (!isPTest || qdMax[snp->GetLabelCount() - 1].Evaluate(sumD)) 
#endif
		{
			//Risk Cells are built, let's build the score for these pieces
			for (uint j=0; j<foldCount; j++) {
				//Let's skip the one that is the validation set
				if (i!=j || foldCount == 1) {
					sumD2+=folds[j].CalculateD2(snp);
				}
			}
			validation.results.CalcTrainingT(sumD, sumD2);
			if (!isPTest && Verbose) {
				cout<<i+1<<" "<<setw(10)<<sumD<<setw(10)<<sumD2<<setw(16)<<right<<snp->GetLabel();
				cout<<setw(8)<<setprecision(4)<<validation.results.trainingT<<endl;
			}
			FoldStatistics stats(snp->GetLabel());
			stats.order = snp->GetLabelCount() - 1;
			stats.trainingFitness = validation.results.trainingT;
			results[i].Evaluate(snp, stats);
			overallResults->Append(i, stats);
	
			if (validation.results.trainingT > bestT)
				bestT = validation.results.trainingT;

		}

		validation.ResetCache();
		snp->AppendEvaluation(bestT);

	}

	
	return true;
}


//Evaluate MAtched Odds Ratio based on the parental status array
/**
 * @param snp 
 * @param status 
 * @param highRiskCells[] 
 * @param hrCellCount 
 */
float  EvalBalancedAccuracyPDT::FindMatchedOddsRatio(SnpAligned *snp, FoldType status) {
	//ostream *os = &cout;
	//if (log)
	//	os = log->GetStream();

	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);
	BitSetType highRisk = snp->GetHrCells();

	uint gtCount = snp->CountGenotypes();
 	FoldType lastFold = &status[sibshipStatusCount];
	uint gtID;
	uint affCount = 0;
	uint unaffCount = 0;
	MatchedOddsRatio oddsRatio;
	stringstream ss;

	uint foldID = 0;
	//Walk through the family masks and ONLY the HR cells
	for (FoldType currCC=status; currCC<lastFold; currCC++) {
		uint truePos=0;
		uint trueNeg=0;
		uint falsePos=0;
		uint falseNeg=0;

		for (gtID=1; gtID<gtCount; gtID++) {
			ss<<"fold: "<<foldID<<"\tGenotype: "<<gtID<<"\n";

			curGenotype=snp->GetGenotype(gtID);
			
			affCount = (currCC->affected & curGenotype).count();
			
			unaffCount = (currCC->unaffected & curGenotype). count();
			//cout<<gtID<<"\t"<<curGenotype<<"\n";
			//Check to see if it's a high risk cell
			if (highRisk[gtID]) {
			//if (cellIdx < hrCellCount && gtID == highRiskCells[cellIdx]) {
				//cellIdx++;
				truePos+=affCount;
				falsePos+=unaffCount;
				ss<<"*";
			}
			else {
				falseNeg+=affCount;
				trueNeg+=unaffCount;
			}
			ss<<"\t"<<curGenotype<<"\t "<<curGenotype<<"\n";
			ss<<"\t"<<currCC->affected<<"\t "<<currCC->unaffected<<"\n";
			ss<<"\t"<<(currCC->affected & curGenotype)<<"\t "<<(currCC->unaffected & curGenotype)<<"\n";
			
		}

		if (truePos+falseNeg != falsePos+trueNeg) {
			BitSetType lostData=snp->Verify();
			size_t pos = lostData.find_first();
			cout<<"Unbalanced parental information: "<<truePos<<" "<<trueNeg<<" "<<falsePos<<" "<<falseNeg<<" ";
			while (pos != BitSetType::npos) {
				cout<<pos<<" ";
				pos = lostData.find_next(pos);
			}
			cout<<"\n";
			cout<<ss.str();
			cout.flush();
			return -1.0;
		}

		oddsRatio.AddFamily(truePos, falsePos, falseNeg, trueNeg);
	}
	//oddsRatio.GenerateReport(os);
	return oddsRatio.GetRatio();
							
}
//Evaluate MAtched Odds Ratio based on the parental status array
/**
 * @param snp 
 * @param status 
 * @param highRiskCells[] 
 * @param hrCellCount 
 * @param mask				- This indicates which individuals are to be considered as part of the fold (1s are considered only)
 */
float  EvalBalancedAccuracyPDT::FindMatchedOddsRatio(SnpAligned *snp, FoldType status, CaseControlStatus &mask) {
	//ostream *os = &cout;
	//if (log)
	//	os = log->GetStream();

	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);
	BitSetType highRisk = snp->GetHrCells();

	uint gtCount = snp->CountGenotypes();
 	FoldType lastFold = &status[sibshipStatusCount];
	uint gtID;
	uint affCount = 0;
	uint unaffCount = 0;
	MatchedOddsRatio oddsRatio;
//	stringstream ss;

	uint foldID = 0;
	//Walk through the family masks and ONLY the HR cells
	for (FoldType currCC=status; currCC<lastFold; currCC++) {
		uint truePos=0;
		uint trueNeg=0;
		uint falsePos=0;
		uint falseNeg=0;

		BitSetType aff = currCC->affected & mask.affected;
		BitSetType unaff = currCC->unaffected  & mask.unaffected;
		
		for (gtID=1; gtID<gtCount; gtID++) {
//cout<<".";cout.flush();
//			ss<<"fold: "<<foldID<<"\tGenotype: "<<gtID<<"\n";

			curGenotype=snp->GetGenotype(gtID);
			affCount = (aff & curGenotype).count();
			unaffCount = (unaff & curGenotype). count();

			//cout<<gtID<<"\t"<<curGenotype<<"\n";
			//Check to see if it's a high risk cell
			if (highRisk[gtID]) {
			//if (cellIdx < hrCellCount && gtID == highRiskCells[cellIdx]) {
				//cellIdx++;
				truePos+=affCount;
				falsePos+=unaffCount;
//				ss<<"*";
			}
			else {
				falseNeg+=affCount;
				trueNeg+=unaffCount;
			}
//			ss<<"\t"<<curGenotype<<"\t "<<curGenotype<<"\n";
//			ss<<"\t"<<aff<<"\t "<<unaff<<"\n";
//			ss<<"\t"<<(aff & curGenotype)<<"\t "<<(unaff & curGenotype)<<"\n";
			
		}

		if (truePos+falseNeg != falsePos+trueNeg) {
			BitSetType lostData=snp->Verify();
			uint pos = lostData.find_first();
			cout<<"Unbalanced parental information: "<<truePos<<" "<<trueNeg<<" "<<falsePos<<" "<<falseNeg<<" ";
//			while (pos != BitSetType::npos) {
//				cout<<pos<<" ";
//				pos = lostData.find_next(pos);
//			}
			cout<<"\n";
//			cout<<ss.str();
			cout.flush();
			return -1.0;
		}
		oddsRatio.AddFamily(truePos, falsePos, falseNeg, trueNeg);
	}
	//oddsRatio.GenerateReport(os);
	return oddsRatio.GetRatio();
							
}

ModelStatistics EvalBalancedAccuracyPDT::EvaluateVerbose(SnpAligned *snp) {
	uint totalAffected;
	uint totalUnaffected;
	uint  sumD = 0;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	uint sumD2 =0;
	
	uint gtCount = snp->CountGenotypes();
	BitSetType highRisk(gtCount);
	//float bestT = 0.0;
	ModelStatistics overallStat;

	//Figure out what the overall ratio is
	BitSetType curGenotype;

	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();
	string tag;
	for (uint i=0; i<foldCount; i++) {

		CaseControlStatus training = folds[i].GetOverallStatus();
		CaseControlStatus testing = folds[i].GetOverallStatus();

		sumD2 = 0;
		sumD  = 0;

		PdtFold &validation = folds[i];
		
		if (foldCount > 1)
			training = validation.GetTrainingSet();
			
		curGenotype = snp->GetGenotype(0);






		//cout<<"Training Mask: "<<training.total<<"\n";
		//cout<<"Testing Mask:  "<<testing.total<<"\n";

		cout<<"Total Individuals in Training Set: "<<training.total.count()<<"\n";
		cout<<"Total Individuals in Testing Set:  "<<testing.total.count()<<"\n";      


		//Figure out what the overall ratio is
		BitSetType curGenotype = snp->GetGenotype(0);
	
		totalAffected = (curGenotype & training.affected).count();
		totalUnaffected = (curGenotype & training.unaffected).count();
		*os<<"\n\nModel Details\n\tModel ID: "<<snp->GetLabel();

		*os<<" Fold "<<i;

		*os<<"\n";
		float totalIndividuals = (float)training.affected.count();

		*os<<"\tTotal Affected: "<<(int)totalAffected;
		*os<<" ("<<setw(6)<<setiosflags(ios::fixed | ios::showpoint);
		*os<<setprecision(2) <<((float)totalAffected/totalIndividuals)*100.0 <<"%)";
		*os<<"\t Total Unaffected: "<<(int)totalUnaffected;
	
		*os<<" ("<<setw(6)<<setiosflags(ios::fixed | ios::showpoint); 
		*os<<setprecision(2)<<((float)totalUnaffected/totalIndividuals)*100.0<<"%)"<<endl;
				
		*os<<"\n          Genotype                 Individual Counts\n";
		stringstream gtLabel;
		for (uint n=0; n<snp->GetLabelCount(); n++) 
			gtLabel<<setw(4)<<right<<snp->GetLabel(n)<<" ";

		*os<<setw(18)<<right<<gtLabel.str();
		*os<<"  "<<setw(12)<<"Affected"<<setw(12)<<"Unaffected";
		*os<<setw(12)<<"Total";
		*os<<setw(12)<<"Ratio";
		*os<<"\n   --------------------------------------------------------------------\n";            

		highRisk.reset();

		for (uint gtID=1; gtID<gtCount; gtID++) {
	
			curGenotype=snp->GetGenotype(gtID);
			//Determine if the cells are high risk
			localAffected = (curGenotype & training.affected).count();
			localUnaffected = (curGenotype & training.unaffected).count();
	
			if (localUnaffected >0) {
				localRatio = (float)localAffected/(float)localUnaffected;		
			}
			else
				localRatio=localAffected;
			
			if (localRatio >= 1.0) {
				tag=" *";
				sumD+=(localAffected - localUnaffected);
				highRisk[gtID]= true;
			}		
			else {
				tag="  ";
			}
			*os<<setw(18)<<right<<snp->GetGenotypeLabel(gtID)<<tag<<setw(12)<<(uint)localAffected;
			*os<<setw(12)<<(uint)localUnaffected<<setw(12)<<((uint)localAffected+(uint)localUnaffected);
			*os<<setw(12)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(4)<<localRatio<<endl;
		}
		snp->SetHrCells( highRisk );
		*os<<"\t\t*Indicates models that are determined to be High Risk\n";
		*os<<"\n\t\t Heterozygote alleles are not necessarily ordered the\n";
		*os<<"\t\t same as they were found in the original data\n";

		//Risk Cells are built, let's build the score for these pieces
		for (uint j=0; j<foldCount; j++) {
			//Let's skip the one that is the validation set
			if (i!=j || foldCount == 1) {
				sumD2+=folds[j].CalculateD2(snp);
			}
		}

		validation.results.CalcTrainingT(sumD, sumD2);
		validation.CalculateAccuracy(snp);

		float matchedOddsRatio = 0.0;
		float classError = validation.results.ClassificationError(); 
		float predError = validation.results.PredictionError();

		*os<<"\n\nClassification Details:\n";
		*os<<setw(12)<<""<<setw(18)<<" Correctly"<<setw(18)<<"Incorrectly"<<endl;
		*os<<setw(12)<<""<<setw(18)<<"Classified"<<setw(18)<<"Classified"<<endl;
		*os<<setw(12)<<"Affected"<<setw(6)<<validation.results.classification.truePos<<" ('High-Risk')";
		*os<<setw(6)<<validation.results.classification.falseNeg<<" ('Low-Risk') "<<endl;
		*os<<setw(12)<<"Unaffected"<<setw(6)<<validation.results.classification.trueNeg<<" ('Low-Risk') ";
		*os<<setw(6)<<validation.results.classification.falsePos<<" ('High-Risk')"<<endl;
 		*os<<setw(14)<<"   "<<setprecision(2)<<100.0*(float)(validation.results.classification.truePos+validation.results.classification.trueNeg)/((float)(validation.results.classification.total) )<<"%";
		*os<<setw(14)<<"   "<<setprecision(2)<<classError<<"%";
		*os<<setw(14)<<" of "<<validation.results.classification.total<<endl;
		*os<<"\nSummary: \nModel ID:"<<snp->GetLabel()<<endl;
		*os<<setw(25)<<right<<"Sum(D): "<<sumD<<endl;
		*os<<setw(25)<<right<<"Sum(D*D): "<<sumD2<<endl;

		*os<<setw(25)<<setprecision(4)<<right<<"T-Statistic: "<<validation.results.trainingT<<endl;

		if (foldCount > 1) 
			*os<<setw(25)<<right<<"Pred. Error: "<<predError<<endl;
		*os<<setw(25)<<right<<"Matched Odds Ratio: "<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3);

		if (sibshipStatus.size() > 0) 
			matchedOddsRatio = FindMatchedOddsRatio(snp, sibshipStatus[0], training);
		*os<<left<<matchedOddsRatio<<endl;

		FoldStatistics fs(snp->GetLabel(), snp->GetLabelCount(), validation.results.trainingT, validation.results.testingT, classError, predError, matchedOddsRatio);
		overallStat.Append(i, fs);
		validation.ResetCache();

	}
	*os<<setw(25)<<"Overall Matched Odds Ratio: "<<FindMatchedOddsRatio(snp, sibshipStatus[0]);
	return overallStat;
}





}

