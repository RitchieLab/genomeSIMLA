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


	*os<<setw(width)<<"Model"<<setw(width)<<"T Statstic";
	if (distribution)
		*os<<setw(width)<<"p Value";
	*os<<endl;

	for (uint p=comboStart-1; p<comboEnd; p++){
		for (uint i=0; i<foldCount; i++) {
			if (distribution) {
				*os<<p+1<<"\t";


				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].highTStatistic[p];
				*os<<"     p < "<<distribution->GetPValue(results[i].topModel[p]->GetLabelCount(), results[i].highTStatistic[p])<<"\n";
			}
			else {
				*os<<p+1<<"\t";


				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].highTStatistic[p]<<"\n";
			}
		}
	}
}
void EvalBalancedAccuracyPDT::ReportResults(ostream *os) {
	*os<<"Top Model of the PDT Search:\n";
	uint width = 12;
	*os<<"Model\n";
	*os<<"Size\t";


	*os<<setw(width)<<"Model"<<setw(width)<<"T Statstic";
	if (distribution)
		*os<<setw(width)<<"p Value";
	*os<<endl;
	for (uint p=comboStart-1; p<comboEnd; p++){
		for (uint i=0; i<foldCount; i++) {
			if (distribution) {
				*os<<p+1<<"\t";

 
				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].highTStatistic[p];
				*os<<"     p < "<<distribution->GetPValue(results[i].topModel[p]->GetLabelCount(), results[i].highTStatistic[p])<<"\n";
			}
			else {
				*os<<p+1<<"\t";


				*os<<setw(width)<<results[i].topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i].highTStatistic[p]<<"\n";
			}
		}
	}
}

//Eventually, this will probably just set up the scores
void EvalBalancedAccuracyPDT::PostEvaluation() {
	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();

	stringstream details;
	//Iterate over each of the status objects insde the vector
	for (uint i=0; i<foldCount; i++) {
		if (!isPTest) {


			*os<<setw(20)<<right<<"Model ";
			if (foldCount > 1) 
				*os<<setw(20)<<right<<"PDT T Stat ";
			else
				*os<<setw(10)<<right<<"PDT T";
			*os<<setw(10)<<right<<"Cls.";


			*os<<endl;
			*os<<setw(20)<<right<<"  ID ";
			if (foldCount > 1) {
				*os<<setw(10)<<right<<"Training";
				*os<<setw(10)<<right<<"Testing";
			}
			else {
				*os<<setw(10)<<right<<"Stat";
			}
			*os<<setw(10)<<right<<"Error";


			*os<<endl;
			*os<<"    ------------------------------------------------\n";
		}
		for (uint p=comboStart-1; p<comboEnd; p++){
			if (isPTest) {
				//cout<<results[i].highTStatistic.size();
				distribution->Append(results[i].highTStatistic[p], testNumber, results[i].topModel[p]);	
			}
			else {
				EvaluateVerbose(i, results[i].topModel[p], details);
			}
		}
		
	}
	
	if (!isPTest) {
		details<<"\t\t*Indicates models that are determined to be High Risk\n";
		details<<"\n\t\t Heterozygote alleles are not necessarily ordered the\n";
		details<<"\t\t same as they were found in the original data\n";
		*os<<details.str()<<"\n";

	}
}


float EvalBalancedAccuracyPDT::EvaluateVerbose(uint fold, SnpAligned *snp, stringstream& details) {
	uint totalAffected;
	uint totalUnaffected;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	uint cellCount=0;
	float bestT = 0.0;
	int  sumD = 0;
	int sumD2 =0;
	int localD = 0;						//This is the D that represents the local slice only
	int localD2 = 0;						//This represents the D2 for the local slice
	
	uint gtCount = snp->CountGenotypes();
	uint highRiskCells[gtCount];

	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();
	string tag;
	sumD2 = 0;
	sumD  = 0;
	cellCount = 0;

	PdtFold &validation = folds[fold];
	CaseControlStatus training = validation.GetOverallStatus();
	CaseControlStatus testing  = validation.GetOverallStatus();
	
	if (foldCount > 1)
		training = validation.GetTrainingSet();
		
	//Figure out what the overall ratio is
	BitSetType curGenotype = ~snp->GetGenotype(0);

	totalAffected = (curGenotype & training.affected).count();
	totalUnaffected = (curGenotype & training.unaffected).count();

	float totalIndividuals = (float)training.affected.size();
	details<<"\n\nModel Details [ "<<snp->GetLabel()<<" ]\n";
	details<<"\t    Affected: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<((float)totalAffected/totalIndividuals)*100.0 <<"%\t";
	details<<"\tUnaffected: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<((float)totalUnaffected/totalIndividuals)*100.0<<"%\n";
	details<<"\tMissing Data: "<<right<<setw(6)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<100.0 - ((float)curGenotype.count()/(float)curGenotype.size())*100.0<<"%\n";
			
	details<<setw(20)<<right<<"\nGenotype"<<setw(18)<<right<<"Individuals      ";
	stringstream gtLabel;
	for (uint n=0; n<snp->GetLabelCount(); n++) 
		gtLabel<<setw(4)<<right<<snp->GetLabel(n)<<" ";

	details<<"\n"<<setw(20)<<right<<gtLabel.str();
	details<<"  "<<setw(6)<<"A"<<setw(6)<<"U";
	details<<setw(6)<<"TOT"<<setw(8)<<"Ratio"<<endl;            
	details<<"    ------------------------------------------------\n";

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
			localD+=((curGenotype & testing.affected).count() - (curGenotype & testing.unaffected).count());
			highRiskCells[cellCount++]= gtID;
		}		
		else {
			tag="  ";
		}
		details<<setw(20)<<right<<snp->GetGenotypeLabel(gtID)<<tag<<setw(6)<<(uint)localAffected;
		details<<setw(6)<<(uint)localUnaffected<<setw(6)<<((uint)localAffected+(uint)localUnaffected);
		details<<setw(8)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(2)<<localRatio<<endl;
	}

	//Risk Cells are built, let's build the score for these pieces
	for (uint j=0; j<foldCount; j++) {
		//Let's skip the one that is the validation set
		if (fold != j || foldCount == 1) {
			sumD2+=folds[j].CalculateD2(highRiskCells, cellCount, snp);
		}
	}

	localD2 = folds[fold].CalculateD2(highRiskCells, cellCount, snp);

	validation.results.CalcTrainingT(sumD, sumD2);
	validation.CalculateAccuracy(highRiskCells, cellCount, snp);
	validation.results.CalcTestingT(localD, localD2);
	
	*os<<setw(19)<<right<<snp->GetLabel()<<" ";
	*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<right<<validation.results.trainingT;
	if (foldCount > 1)
		*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<right<<validation.results.testingT;
	
	*os<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<right;
	*os<<100.0*(float)(validation.results.falsePos+validation.results.falseNeg)/((float)validation.results.total)<<"%";


	*os<<endl;
	if (validation.results.trainingT > bestT)
		bestT = validation.results.trainingT;

	return bestT;
}

inline
bool EvalBalancedAccuracyPDT::EvaluateModel(SnpAligned *snp, uint position) {
	uint  sumD = 0;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	uint sumD2 =0;
	
	uint gtCount = snp->CountGenotypes();
	uint highRiskCells[gtCount];
	uint cellCount=0;
	float bestT = 0.0;

	//Figure out what the overall ratio is
	BitSetType curGenotype;

	CaseControlStatus trainingSet = folds[0].GetOverallStatus();
	for (uint i=0; i<foldCount; i++) {
		sumD2 = 0;
		sumD  = 0;
		cellCount = 0;

		PdtFold &validation = folds[i];
		if (foldCount > 1) 
			trainingSet = validation.GetTrainingSet();
		curGenotype = snp->GetGenotype(0);
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
				highRiskCells[cellCount++]= gtID;
			}
				
		}
		//Risk Cells are built, let's build the score for these pieces
		for (uint j=0; j<foldCount; j++) {
			//Let's skip the one that is the validation set
			if (i!=j || foldCount == 1) {
				sumD2+=folds[j].CalculateD2(highRiskCells, cellCount, snp);
			}
		}
		validation.results.CalcTrainingT(sumD, sumD2);
		if (!isPTest && Verbose) {
			cout<<i+1<<" "<<setw(10)<<sumD<<setw(10)<<sumD2<<setw(16)<<right<<snp->GetLabel();
			cout<<setw(8)<<setprecision(4)<<validation.results.trainingT<<endl;
		}
		results[i].Evaluate(snp, validation.results);

		if (validation.results.trainingT > bestT)
			bestT = validation.results.trainingT;
	}

	snp->SetLastEvaluation(bestT);
	
	return true;
}


//Evaluate MAtched Odds Ratio based on the parental status array
/**
 * @param snp 
 * @param status 
 * @param highRiskCells[] 
 * @param hrCellCount 
 */
float  EvalBalancedAccuracyPDT::FindMatchedOddsRatio(SnpAligned *snp, FoldType status, uint highRiskCells[], uint hrCellCount) {
	//ostream *os = &cout;
	//if (log)
	//	os = log->GetStream();

	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);
	uint gtCount = snp->CountGenotypes();
 	FoldType lastFold = &status[sibshipStatusCount];
	uint gtID;
	uint affCount = 0;
	uint unaffCount = 0;
	MatchedOddsRatio oddsRatio;
	//Walk through the family masks and ONLY the HR cells
	for (FoldType currCC=status; currCC<lastFold; currCC++) {
		uint truePos=0;
		uint trueNeg=0;
		uint falsePos=0;
		uint falseNeg=0;

		uint cellIdx = 0;

		for (gtID=1; gtID<gtCount; gtID++) {
			curGenotype=snp->GetGenotype(gtID);
			affCount = (currCC->affected & curGenotype).count();
			unaffCount = (currCC->unaffected & curGenotype). count();
			//cout<<gtID<<"\t"<<curGenotype<<"\n";
			//Check to see if it's a high risk cell
			if (cellIdx < hrCellCount && gtID == highRiskCells[cellIdx]) {
				cellIdx++;
				truePos+=affCount;
				falsePos+=unaffCount;
			}
			else {
				falseNeg+=affCount;
				trueNeg+=unaffCount;
			}
			
		}

		if (truePos+falseNeg != falsePos+trueNeg) {
			BitSetType lostData=snp->Verify();
			uint pos = lostData.find_first();
			//cout<<"Unbalanced parental information: "<<truePos<<" "<<trueNeg<<" "<<falsePos<<" "<<falseNeg<<" ";
			while (pos != BitSetType::npos) {
				cout<<pos<<" ";
				pos = lostData.find_next(pos);
			}
			//cout<<"\n";
			return -1.0;
		}

		oddsRatio.AddFamily(truePos, falsePos, falseNeg, trueNeg);
	}
	//oddsRatio.GenerateReport(os);
	return oddsRatio.GetRatio();
							
}


float EvalBalancedAccuracyPDT::EvaluateVerbose(SnpAligned *snp) {
	uint totalAffected;
	uint totalUnaffected;
	uint  sumD = 0;

	uint localAffected;
	uint localUnaffected;
	float localRatio;
	uint sumD2 =0;
	
	uint gtCount = snp->CountGenotypes();
	uint highRiskCells[gtCount];
	uint cellCount=0;
	float bestT = 0.0;

	//Figure out what the overall ratio is
	BitSetType curGenotype;

	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();
	string tag;
	CaseControlStatus training = folds[0].GetOverallStatus();
	for (uint i=0; i<foldCount; i++) {
		sumD2 = 0;
		sumD  = 0;
		cellCount = 0;

		PdtFold &validation = folds[i];
		
		if (foldCount > 1)
			training = validation.GetTrainingSet();
			
		curGenotype = snp->GetGenotype(0);

		//Figure out what the overall ratio is
		BitSetType curGenotype = ~snp->GetGenotype(0);
	
		totalAffected = (curGenotype & training.affected).count();
		totalUnaffected = (curGenotype & training.unaffected).count();
		*os<<"\n\nModel Details\n\tModel ID: "<<snp->GetLabel();


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
				highRiskCells[cellCount++]= gtID;
			}		
			else {
				tag="  ";
			}
			*os<<setw(18)<<right<<snp->GetGenotypeLabel(gtID)<<tag<<setw(12)<<(uint)localAffected;
			*os<<setw(12)<<(uint)localUnaffected<<setw(12)<<((uint)localAffected+(uint)localUnaffected);
			*os<<setw(12)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(4)<<localRatio<<endl;
		}
		*os<<"\t\t*Indicates models that are determined to be High Risk\n";
		*os<<"\n\t\t Heterozygote alleles are not necessarily ordered the\n";
		*os<<"\t\t same as they were found in the original data\n";

		//Risk Cells are built, let's build the score for these pieces
		for (uint j=0; j<foldCount; j++) {
			//Let's skip the one that is the validation set
			if (i!=j || foldCount == 1) {
				sumD2+=folds[j].CalculateD2(highRiskCells, cellCount, snp);
			}
		}

		validation.results.CalcTrainingT(sumD, sumD2);
		validation.CalculateAccuracy(highRiskCells, cellCount, snp);

		*os<<"\n\nClassification Details:\n";
		*os<<setw(12)<<""<<setw(18)<<" Correctly"<<setw(18)<<"Incorrectly"<<endl;
		*os<<setw(12)<<""<<setw(18)<<"Classified"<<setw(18)<<"Classified"<<endl;
		*os<<setw(12)<<"Affected"<<setw(6)<<validation.results.truePos<<" ('High-Risk')";
		*os<<setw(6)<<validation.results.falseNeg<<" ('Low-Risk') "<<endl;
		*os<<setw(12)<<"Unaffected"<<setw(6)<<validation.results.trueNeg<<" ('Low-Risk') ";
		*os<<setw(6)<<validation.results.falsePos<<" ('High-Risk')"<<endl;
 		*os<<setw(14)<<"   "<<setprecision(2)<<100.0*(float)(validation.results.truePos+validation.results.trueNeg)/((float)(validation.results.total) )<<"%";
		*os<<setw(14)<<"   "<<setprecision(2)<<100.0*(float)(validation.results.falsePos+validation.results.falseNeg)/((float)validation.results.total)<<"%";
		*os<<setw(14)<<" of "<<validation.results.total<<endl;
		*os<<"\nSummary: \nModel ID:"<<snp->GetLabel()<<endl;
		*os<<setw(25)<<right<<"Sum(D): "<<sumD<<endl;
		*os<<setw(25)<<right<<"Sum(D*D): "<<sumD2<<"\n";
		*os<<setw(25)<<setprecision(4)<<right<<"T-Statistic: "<<validation.results.trainingT<<endl;

		
		if (validation.results.trainingT > bestT)
			bestT = validation.results.trainingT;

	}
	return bestT;
}





}
