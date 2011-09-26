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
namespace ESE {

bool EvalBalancedAccuracyPDT::Verbose = false;

void EvalBalancedAccuracyPDT::SummarizeResults(ostream *os) {
	uint width = 12;
	*os<<setw(width)<<"Model"<<setw(width)<<"T Statstic";
	if (distribution)
		*os<<setw(width)<<"p Value";
	*os<<endl;
	uint count=results.size();
	for (uint i=0; i<count; i++) {
		for (uint p=comboStart-1; p<comboEnd; p++){
			if (distribution) {
				*os<<setw(width)<<results[i]->topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i]->highTStatistic[p];
				*os<<"     p < "<<distribution->GetPValue(results[i]->topModel[p]->GetLabelCount(), results[i]->highTStatistic[p])<<"\n";
			}
			else {
				*os<<setw(width)<<results[i]->topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i]->highTStatistic[p]<<"\n";
			}
		}
	}
}
void EvalBalancedAccuracyPDT::ReportResults(ostream *os) {
	*os<<"Top Model of the PDT Search:\n";
	uint width = 12;
	*os<<setw(width)<<"Model"<<setw(width)<<"T Statstic";
	if (distribution)
		*os<<setw(width)<<"p Value";
	*os<<endl;
	uint count=results.size();
	for (uint i=0; i<count; i++) {
		for (uint p=comboStart-1; p<comboEnd; p++){
			if (distribution) {
				*os<<setw(width)<<results[i]->topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i]->highTStatistic[p];
				*os<<"     p < "<<distribution->GetPValue(results[i]->topModel[p]->GetLabelCount(), results[i]->highTStatistic[p])<<"\n";
			}
			else {
				*os<<setw(width)<<results[i]->topModel[p]->GetLabel();
				*os<<setw(width)<<setprecision(4)<<results[i]->highTStatistic[p]<<"\n";
			}
		}
	}
}

//Eventually, this will probably just set up the scores
void EvalBalancedAccuracyPDT::PostEvaluation() {
	uint count=results.size();
	//Iterate over each of the status objects insde the vector
	for (uint i=0; i<count; i++) {
		for (uint p=comboStart-1; p<comboEnd; p++){
			FoldType ccSet = testingSet[i];
			if (isPTest) {
				Results *res=results[i];
				if (res) {
					cout<<res->highTStatistic.size();
					distribution->Append(res->highTStatistic[p], testNumber, res->topModel[p]);	
				}
				else {
					cout<<"Trying to read into a result, "<<i<<", object that doesn't exist\n";
					assert(0);
				}
			}
			else
				EvaluateVerbose(results[i]->topModel[p], ccSet);
		}
	}
}

//Evaluate MAtched Odds Ratio based on the parental status array
/**
 * @param snp 
 * @param status 
 * @param highRiskCells[] 
 * @param hrCellCount 
 */
void EvalBalancedAccuracyPDT::FindMatchedOddsRatio(SnpAligned *snp, FoldType status, uint highRiskCells[], uint hrCellCount) {
	ostream *os = &cout;
	if (log)
		os = log->GetStream();

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
				falseNeg+=unaffCount;
			}
			else {
				falsePos+=affCount;
				trueNeg+=unaffCount;
			}
			
		}

		if (truePos+falsePos != falseNeg+trueNeg) {
			BitSetType lostData=snp->Verify();
			uint pos = lostData.find_first();
			cout<<"Unbalanced parental information: ";
			while (pos != BitSetType::npos) {
				cout<<pos<<" ";
				pos = lostData.find_next(pos);
			}
			cout<<"\n";
			return;
		}

		oddsRatio.AddFamily(truePos, falsePos, falseNeg, trueNeg);
	}
	oddsRatio.GenerateReport(os);
							
}

float EvalBalancedAccuracyPDT::EvaluateVerbose(SnpAligned *snp, FoldType fold) {
	ostream *os = &cout;
	if (log && log->GetStream())
		os = log->GetStream();

	uint gtCount = snp->CountGenotypes();

	//Initialize the vector to be the size of the
	float localD = 0.0; 
	float sumD = 0.0;
	float sumD2 = 0.0;

	//Let's ID the genotypes that are HR and then only scan those
	float affCount = 0.0;
	float unaffCount = 0.0;
	float totalAffected = 0.0;
	float totalUnaffected = 0.0;
	float localAffected = 0.0;
	float localUnaffected = 0.0;
	float totalRatio= 0.0;
	float localRatio =0.0;
	float tStatistic=0.0;
	uint grandTotal =0;
	uint truePos=0;
	uint trueNeg=0;
	uint falsePos=0;
	uint falseNeg=0;
	string tag;
	uint highRiskCells[gtCount];
	uint cellCount=0;

	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);
	
	totalAffected = (float)(curGenotype & overallStatus.affected).count();
	totalUnaffected = (float)(curGenotype & overallStatus.unaffected).count();
	grandTotal = (uint)(totalAffected + totalUnaffected);
	*os<<"\n\nModel Details\n\tModel ID: "<<snp->GetLabel();

	*os<<setw(30)<<"\n\tTotal Affected: "<<(int)totalAffected;
	*os<<"("<<setw(6)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3) <<((float)totalAffected / (float)overallStatus.affected.count())*100.0 <<"%)";
	*os<<setw(25)<<"  Total Unaffected: "<<(int)totalUnaffected;
	
	*os<<"("<<setw(6)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3) <<((float)totalUnaffected/(float)overallStatus.unaffected.count())*100.0<<"%)"<<endl;
	totalRatio = totalAffected  / totalUnaffected;

	if (!useModelThreshold) {
		totalRatio = overallRatio;
		*os<<setw(30)<<"Application Threshold: "<<totalRatio<<endl;
	} 

	*os<<"          Genotype ";

	//Set up our model details
	stringstream gtLabel;
	
	for (uint i=0; i<snp->GetLabelCount(); i++) 
		gtLabel<<setw(4)<<right<<snp->GetLabel(i)<<" ";

	*os<<"\n"<<setw(18)<<right<<gtLabel.str()<<"  "<<setw(12)<<"Affected"<<setw(12)<<"Unaffected"<<setw(12)<<"Total"<<setw(12)<<"Ratio"<<endl;            
	for (uint gtID=1; gtID<gtCount; gtID++) {
		curGenotype=snp->GetGenotype(gtID);
	
		//Determine if the cells are high risk
		localAffected = (float)(curGenotype & overallStatus.affected).count();
		localUnaffected = (float)(curGenotype & overallStatus.unaffected).count();

		if (localUnaffected > 0.0) {
			localRatio = localAffected/localUnaffected;		
		}
		else
			localRatio=localAffected;
		
		if (localRatio >= totalRatio) {
			tag=" *";
			sumD+=(localAffected - localUnaffected);
			highRiskCells[cellCount++]= gtID;
			truePos+=(uint)localAffected;
			falseNeg+=(uint)localUnaffected;
		}
		else {
			tag="  ";
			falsePos+=(uint)localAffected;
			trueNeg+=(uint)localUnaffected;
		}
		*os<<setw(18)<<right<<snp->GetGenotypeLabel(gtID)<<tag<<setw(12)<<(uint)localAffected;
		*os<<setw(12)<<(uint)localUnaffected<<setw(12)<<((uint)localAffected+(uint)localUnaffected);
		*os<<setw(12)<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(4)<<localRatio<<endl;
	}
	*os<<"\t\t*Indicates models that are determined to be High Risk\n";
	*os<<"\n\t\t Heterozygote alleles are not necessarily ordered the\n";
	*os<<"\t\t same as they were found in the original data\n";

	*os<<"\n\nClassification Details:\n";
	*os<<setw(12)<<""<<setw(18)<<" Correctly"<<setw(18)<<"Incorrectly"<<endl;
	*os<<setw(12)<<""<<setw(18)<<"Classified"<<setw(18)<<"Classified"<<endl;
	*os<<setw(12)<<"Affected"<<setw(6)<<truePos  <<" ('High-Risk')"<<setw(6)<<falsePos<<" ('Low-Risk') "<<endl;
	*os<<setw(12)<<"Unaffected"<<setw(6)<<trueNeg<<" ('Low-Risk') "<<setw(6)<<falseNeg<<" ('High-Risk')"<<endl;
 	*os<<setw(14)<<" "<<setprecision(4)<<100.0*(float)(truePos+trueNeg)/((float)(grandTotal) )<<"%";
	*os<<setw(12)<<" "<<setprecision(4)<<100.0*(float)(falsePos+falseNeg)/((float)grandTotal)<<"%"<<setw(12)<<" of "<<grandTotal<<endl;

	
	FoldType lastFold = &fold[statusCount];
	uint gtID;
	
	//This isn't ideal

	FindMatchedOddsRatio(snp, sibshipStatus[0], highRiskCells, cellCount);
	//*os<<"\nD={ ";

	truePos=0;
	trueNeg=0;
	falsePos=0;
	falseNeg=0;
	//Walk through the family masks and ONLY the HR cells
	for (FoldType currCC=fold; currCC<lastFold; currCC++) {

		localD = 0;
		uint cellIdx = 0;
		for (gtID=1; gtID<gtCount; gtID++) {
			curGenotype=snp->GetGenotype(gtID);
			
			affCount = (float)(currCC->affected & curGenotype).count();
			unaffCount = (float)(currCC->unaffected & curGenotype). count();

			//Check to see if it's a high risk cell
			if (cellIdx < cellCount && gtID == highRiskCells[cellIdx]) {
				cellIdx++;
				localD+=(affCount - unaffCount);
			}
		}
		//Sum up Ds for each family
		//*os<<(int)localD<<" ";
		sumD2+= localD * localD;

	}					
	//*os<<" }\n";
	
	if (sumD2 > 0)
		tStatistic=sumD/(sqrt(sumD2));
	*os<<"Summary: \n  Model ID: "<<snp->GetLabel()<<endl;
	*os<<setw(25)<<right<<"Sum(D): "<<(uint)sumD<<endl;
	*os<<setw(25)<<right<<"Sum(D*D): "<<(uint)sumD2<<"\n";
	*os<<setw(25)<<right<<"T-Statistic: "<<tStatistic<<endl;

	snp->SetLastEvaluation(tStatistic);
	return tStatistic;
}



float EvalBalancedAccuracyPDT::Evaluate(SnpAligned *snp, FoldType fold) {
	//uint statusCount = fold.size();
	uint gtCount = snp->CountGenotypes();

	//Initialize the vector to be the size of the
	float localD = 0.0; 
	float sumD = 0.0;
	float otherSumD = 0.0;
	float sumD2 = 0.0;

	//Let's ID the genotypes that are HR and then only scan those

	float affCount = 0.0;
	float unaffCount = 0.0;
	float totalAffected = 0.0;
	float totalUnaffected = 0.0;
	float localAffected = 0.0;
	float localUnaffected = 0.0;
	float totalRatio= 0.0;
	float localRatio =0.0;
	float tStatistic=0.0;

	uint highRiskCells[gtCount];
	uint cellCount=0;
	//vector<uint> highRiskCells;

	//Figure out what the overall ratio is
	BitSetType curGenotype = snp->GetGenotype(0);

	if (!useModelThreshold)  
		totalRatio = overallRatio;
	else {
		totalAffected = (float)(curGenotype & overallStatus.affected).count();
		totalUnaffected = (float)(curGenotype & overallStatus.unaffected).count();

		//If our data is so far
		if (totalAffected * totalUnaffected == 0) {
			if (!isPTest)
				cout<<setw(16)<<right<<snp->GetLabel()<<setw(8)<<0<<setw(8)<<setprecision(3)<<0<<endl;
			snp->SetLastEvaluation(0.0);
			return 0.0;
		}
		totalRatio = totalAffected  / totalUnaffected;
	}
	
	for (uint gtID=1; gtID<gtCount; gtID++) {
		curGenotype=snp->GetGenotype(gtID);
		//Determine if the cells are high risk
		localAffected = (float)(curGenotype & overallStatus.affected).count();
		localUnaffected = (float)(curGenotype & overallStatus.unaffected).count();

		if (localUnaffected > 0.0) {
			localRatio = localAffected/localUnaffected;		
		}
		else
			localRatio=localAffected;
		
		if (localRatio >= totalRatio) {
			sumD+=(localAffected - localUnaffected);
			highRiskCells[cellCount++]= gtID;
		}
			
	}
	FoldType lastFold = &fold[statusCount];
#ifdef USE_DOPTIMIZATION
	if (qdMax[snp->GetLabelCount()-1].Evaluate(sumD) ){
#endif	
		uint gtID;
		//Walk through the family masks
		for (FoldType currCC=fold; currCC<lastFold; currCC++) {
			//affected=currCC->affected;
			//unaffected=currCC->unaffected;
	
			localD = 0;
			for (gtID=0; gtID<cellCount; gtID++) {
				curGenotype=snp->GetGenotype(highRiskCells[gtID]);
				
				affCount = (float)(currCC->affected & curGenotype).count();
				unaffCount = (float)(currCC->unaffected & curGenotype). count();
				localD+=(affCount - unaffCount);
	
			}
			//Sum up Ds for each family
			otherSumD += localD;
			sumD2+= localD * localD;
		}					
		
		if (sumD > 0)
			tStatistic=sumD/(sqrt(sumD2));

		if (!isPTest && Verbose)
			cout<<setw(16)<<right<<snp->GetLabel()<<setw(8)<<(uint)sumD<<setw(8)<<(uint)sumD2<<setw(8)<<(uint)otherSumD<<setw(8)<<setprecision(4)<<tStatistic<<endl;
#ifdef USE_DOPTIMIZATION
	}
	else
		if (!isPTest && Verbose)
			cout<<setw(16)<<right<<snp->GetLabel()<<setw(8)<<sumD<<"  - Optimized out\n";
#endif
	snp->SetLastEvaluation(tStatistic);
	return tStatistic;
}



}
