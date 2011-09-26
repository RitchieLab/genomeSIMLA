//
// C++ Implementation: configurationparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "evalbalancedaccuracy.h"
#include "evalbalancedaccuracypdt.h"
#include "evalmaxdifference.h"
#include "eseconfiguration.h"
#include "genetics/ptestdistribution.h"
#include "genetics/familynode.h"
#include "matchedoddsratio.h"
#include "snpsearchapplication.h"









//using namespace Genetics::Parser;
using namespace Genetics::Reporting;
using namespace Genetics::Evaluation;

namespace MDR {

////Mendelian error catching/stripping settings
char *EseConfiguration::ConfigFileValues::MendelAction 				= "MENDELIAN_ERROR_LEVEL";
char *EseConfiguration::ConfigFileValues::MendelLocThresh			= "MENDELIAN_LOCUS_THRESHOLD";
char *EseConfiguration::ConfigFileValues::MendelPedThresh			= "MENDELIAN_PEDIGREE_THRESHOLD";
char *EseConfiguration::ConfigFileValues::WriteCleanDatafile		= "WRITE_CLEAN_DATAFILE";


////Fundamental Settings
									//The first locus count for models of interest
char *EseConfiguration::ConfigFileValues::ComboStart				= "COMBO_START";
									//The largest size of interest
char *EseConfiguration::ConfigFileValues::ComboEnd					= "COMBO_END";			
char *EseConfiguration::ConfigFileValues::CalcThreshold				= "USE_MODEL_THRESHOLD";
char *EseConfiguration::ConfigFileValues::AffectedValue 			= "AFFECTED_VALUE";
char *EseConfiguration::ConfigFileValues::UnaffectedValue			= "UNAFFECTED_VALUE";
char *EseConfiguration::ConfigFileValues::LociExclusionList 		= "EXCLUDE_LOCUS";

////Settings for overrideing various report extensions
									//Let user set the extension for distribution reports
char *EseConfiguration::ConfigFileValues::ExtDistribution			= "EXT_DISTRIBUTION";
									//Let user set the extension for the overall report
char *EseConfiguration::ConfigFileValues::ExtReport					= "EXT_REPORT";
									//Let user set the extension used for locus reports		
char *EseConfiguration::ConfigFileValues::ExtLocus					= "EXT_LOCUS";
									//Let user set the extension used for the pedigree report
char *EseConfiguration::ConfigFileValues::ExtPedigree				= "EXT_PEDIGREE";
char *EseConfiguration::ConfigFileValues::ExtCleanedData 			= "EXT_CLEANED_DATA";

									//Let user set the report name
char *EseConfiguration::ConfigFileValues::ReportName				= "REPORT_NAME";					

////Overriding other reporting details
								//Cutoff for the PValue output
char *EseConfiguration::ConfigFileValues::MinPVal					= "PVAL_THRESHOLD";		
								//Number of PTests
char *EseConfiguration::ConfigFileValues::PTestCount				= "PTEST_COUNT";
								//The seed used for the runs
char *EseConfiguration::ConfigFileValues::PTestSeed					= "PTEST_SEED";
								//Used to set the number of cross validation folds
char *EseConfiguration::ConfigFileValues::CrossValFoldCount			= "CROSSVALINTERVAL";
	
//Overriding input details
									//<Used to override the default input format
char *EseConfiguration::ConfigFileValues::_InputFile				="INPUTFILE";
									//Used to override the default format of input data
char *EseConfiguration::ConfigFileValues::_InputFormat				="INPUTFORMAT";
									//<Used to specify to use MDR formatted input
char *EseConfiguration::ConfigFileValues::_MdrFormat				="MDRFORMAT";
									//Used to specify pedigree formatted input
char *EseConfiguration::ConfigFileValues::PedigreeFormat			="PEDIGREE";

//Overriding how the data is to be interpretted
									//<Use to override the max genotypes
char *EseConfiguration::ConfigFileValues::_GenotypeCount			="GENOTYPES";
									//<Used to override the type of analysis
char *EseConfiguration::ConfigFileValues::_AnalysisStyle			="ANALYSISSTYLE";
									//<Used to specify Max Difference on CC data			
char *EseConfiguration::ConfigFileValues::_MaxDiff					="MAXIMIZEDIFF";
									//<Used to specify Balanced Accuracy Analyses on CC data
char *EseConfiguration::ConfigFileValues::BalancedAccuracy 			="BALANCEDACCURACY";
									//<Used to specify PDT analyses on pedigree data
char *EseConfiguration::ConfigFileValues::PDTAnalysis				="PDT";
		
//MD Specific Settings
									//Set the number of models reported on
char *EseConfiguration::ConfigFileValues::_ReportModelCount			="REPORTMODELCOUNT\0";
									//Set the number of affected (required for Base Pair layout)
char *EseConfiguration::ConfigFileValues::_AffectedCount			="AFFECTEDCOUNT\0";
									//Set the number of unaffected (required for Base Pair layout)
char *EseConfiguration::ConfigFileValues::_UnaffectedCount			="UNAFFECTEDCOUNT\0";
									//Set the difference threshold (minimum required to be reported on)
char *EseConfiguration::ConfigFileValues::_BestDifferenceThreshold	="BESTDIFFERENCETHRESHOLD\0";
									//set whether or not the report is sorted
char *EseConfiguration::ConfigFileValues::_SortReportModels			="PERFORMSORT\0";
	
//Pedigree Specific settings
									//<The percentage within the highest D to be considered
char *EseConfiguration::ConfigFileValues::MinD						="DOPTIMIZATION";
char *EseConfiguration::ConfigFileValues::CreateVirtualSibs			="EXPAND_TRIOS";
char *EseConfiguration::ConfigFileValues::ExpandAllAffected			="EXPAND_ALL";
char *EseConfiguration::ConfigFileValues::PedExclusionList			="PEDIGREE_EXCLUSIONS";
char *EseConfiguration::ConfigFileValues::VerboseMOR				="VERBOSE_MATCHED_ODDS_RATIO";
char *EseConfiguration::ConfigFileValues::VerbosePDT				="VERBOSE_EVALUATION";
char *EseConfiguration::ConfigFileValues::WriteFoldDetails			= "VERBOSE_FOLDING";

char *EseConfiguration::ConfigFileValues::CacheD2					="CACHED2";

char *EseConfiguration::ConfigFileValues::WriteRegressionScripts	= "WRITE_REGRESSION_SCRIPTS";


uint EseConfiguration::_instanceCount 								= 0;
EseConfiguration *EseConfiguration::_instance						= NULL;

EseConfiguration::EseConfiguration() : pTestCount(1000), pValueThreshold(0.30), minDThresh(0.0), pTestSeed(13), logLocus(NULL), comboEnd(1),comboStart(1), fDist(NULL), oDist(NULL), eDist(NULL), evalSuite(NULL), inputParser(NULL), logDistribution(NULL), logReport(NULL), logPedigree(NULL), performPerfectSearch(false)
{
	//generateRegressionScript = false;
	families = 	NULL;
	errorMsg="";
	doWriteCleanedDatafile = false;
	cleanedDataset = NULL;

	SetDefaultValues();

	PostLoad();
}


EseConfiguration::~EseConfiguration()
{
	if (inputParser)
		delete inputParser;
	if (logLocus)
		delete logLocus;
	if (evalSuite)
		delete evalSuite;
	if (fDist)
		delete fDist;
	if (oDist)
		delete oDist;
	if (eDist)
		delete eDist;

	//Delete the logs
	if (logDistribution)
		delete logDistribution;
	if (logReport)
		delete logReport;
	if (logPedigree) 
		delete logPedigree;
	if (cleanedDataset)
		delete cleanedDataset;

	if (families)
		FamilyRepository::Release();
}

void EseConfiguration::SetDefaultValues() {
	reportName 				= "";
	extDistribution 		= "pdist";
	extLocus 				= "locus";
	extPedigree 			= "pedigree";
	extCleanedData 			= "ped.cleaned";

	extReport 				= "STDOUT";

	locusWriteStyle			= LocusLogAscii::WriteAll;
	analysisStyle			= PDT;
	reportModelCount		= 0;
	affectedCount			= 0;							//This is illegal- so it will generate an error if it isn't set 
	unaffectedCount 		= 0;
	bestDifferenceThreshold = 0;							//This is illegal if the style is set to one of the max diff styles
	doSortReport			= true;
	inputFile				= "";							//Again, this is illegal 
	inputFormat				= InputSpaceDelimited;
	logLocus				= NULL;							//We don't want to allow this! Output could be a total mess!
	maxGenotypes			= 3;							//We use this to determine who many levels can be expected in the data
	doCreateVirtualSibs		= true;
	doCalcThreshold			= true;
	crossValidationCount 	= 1;
}

void EseConfiguration::PostLoad() {
	if (analysisStyle != PDT && logLocus) {
		delete logLocus;
		logLocus = new LocusLogAscii(locusWriteStyle, reportName.c_str(), extPedigree.c_str());
	}
}

void EseConfiguration::SetInputFormat(const char *value) {
	if (strcmp(value, ConfigFileValues::_MdrFormat)==0) 
		inputFormat=MdrFormat;
	else if (strcmp(value, ConfigFileValues::PedigreeFormat)==0)
		inputFormat=Pedigree;

}

void EseConfiguration::SetAnalysisStyle(const char *value) {
	if (strcmp(value, ConfigFileValues::BalancedAccuracy)==0)
		analysisStyle=BalancedAccuracy;
	else if (strcmp(value, ConfigFileValues::PDTAnalysis)==0)
		analysisStyle=PDT;
	else if (strcmp(value, ConfigFileValues::_MaxDiff)==0) {
		analysisStyle=MaxDiff;
		performPerfectSearch = true;
	}

}



bool EseConfiguration::SetValues(const char *key, const char *value, const char *remainder) {
	bool success=true;
	assert(strlen(key)>0);
	if (strcmp(key, ConfigFileValues::_AnalysisStyle) == 0) 
		SetAnalysisStyle(value);
	else if (strcmp(key, ConfigFileValues::_InputFile) ==0) {
		inputFile=value;
		if (reportName == "") 
			reportName = inputFile;
	}
	else if (strcmp(key, ConfigFileValues::_GenotypeCount)==0)
		maxGenotypes=atoi(value);
	else if (strcmp(key, ConfigFileValues::PTestCount)==0)
		pTestCount=atoi(value);
	else if (strcmp(key, ConfigFileValues::PTestSeed) == 0)
		pTestSeed = atoi(value);
	else if (strcmp(key, ConfigFileValues::MinPVal)==0)
		pValueThreshold=atof(value);
	else if (strcmp(key, ConfigFileValues::MinD)==0) {
		cout<<"Setting Optimization to: "<<value<<" \n";
		minDThresh=atof(value);
	}
	else if (strcmp(key, ConfigFileValues::_InputFormat)==0) 
		SetInputFormat(value);		
	else if (strcmp(key, ConfigFileValues::ComboEnd)==0)
		comboEnd=atoi(value);
	else if (strcmp(key, ConfigFileValues::ComboStart)==0)
		comboStart=atoi(value);
	else if (strcmp(key, ConfigFileValues::_AffectedCount) == 0)
		affectedCount=atoi(value);
	else if (strcmp(key, ConfigFileValues::_UnaffectedCount) == 0)
		unaffectedCount=atoi(value);
	else if (strcmp(key, ConfigFileValues::_BestDifferenceThreshold) == 0)
		bestDifferenceThreshold=atoi(value);
	else if (strcmp(key, ConfigFileValues::_ReportModelCount)==0)
		reportModelCount=atoi(value);
	else if (strcmp(key, ConfigFileValues::_SortReportModels)==0)
		doSortReport=GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::CreateVirtualSibs) == 0)
		FamilyNode::expandTrios = GetBoolean(value);
		//doCreateVirtualSibs=GetBoolean(value);
		
	else if (strcmp(key, ConfigFileValues::ExpandAllAffected) == 0)
		FamilyNode::expandAllAffecteds = GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::AffectedValue) == 0)
		FamilyMember::_affectedValue = atoi(value);
	else if (strcmp(key, ConfigFileValues::UnaffectedValue) == 0)
		FamilyMember::_unaffectedValue = atoi(value);
	else if (strcmp(key, ConfigFileValues::VerboseMOR) ==0)
		MatchedOddsRatio::Verbose=GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::VerbosePDT) == 0)
		EvalBalancedAccuracyPDT::Verbose=GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::CalcThreshold) == 0)
		doCalcThreshold = GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::PedExclusionList) ==0)
		Utility::TokenizeString(remainder, pedExclusions, " ,\t");
	else if (strcmp(key, ConfigFileValues::LociExclusionList) == 0)
		Utility::TokenizeString(remainder, lociExclusions, " ,\t");

	else if (strcmp(key, ConfigFileValues::CrossValFoldCount) == 0)
		crossValidationCount = atoi(value);

	else if (strcmp(key, ConfigFileValues::WriteFoldDetails) ==0) 
		GtLineParserMdrPdt::WriteFoldDetails = GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::ExtReport) == 0) 
		if (strcmp(value, "NOLOG") != 0)
			extReport = value;
		else
			cout<<"* Refusing to execute without returning results. Request to turn analysis results off ignored.\n"; 
	else if (strcmp(key, ConfigFileValues::ExtPedigree) ==0)
		extPedigree = value;
	else if (strcmp(key, ConfigFileValues::ExtDistribution) == 0)
		extDistribution = value;
	else if (strcmp(key, ConfigFileValues::ExtCleanedData) == 0)
		extCleanedData = value;
	else if (strcmp(key, ConfigFileValues::CacheD2) == 0)
		PdtFold::CacheD2 = GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::MendelAction) == 0) 
		mendelianSettings.actionLevel = atoi(value);
	else if (strcmp(key, ConfigFileValues::MendelLocThresh) == 0)
		mendelianSettings.locThreshold = atoi(value);
	else if (strcmp(key, ConfigFileValues::MendelPedThresh) ==0)
		mendelianSettings.pedThreshold = atoi(value);
	else if (strcmp(key, ConfigFileValues::WriteCleanDatafile) == 0)
		doWriteCleanedDatafile = GetBoolean(value);
	else if (strcmp(key, ConfigFileValues::WriteRegressionScripts) == 0)
		SnpSearchApplication::generateRegressionScript = GetBoolean(value);
	else
		success=false;
	return success;
}

ostream *EseConfiguration::GetCleanedDataFile() {
	if (cleanedDataset == NULL)
		cleanedDataset = new AsciiLog(reportName.c_str(), extCleanedData.c_str(), 0);
	return cleanedDataset->GetStream();
}
BasicLog *EseConfiguration::GetDistributionLog() {
	if (logDistribution == NULL) 
		logDistribution = new AsciiLog(reportName.c_str(), extDistribution.c_str(), 0);
	return logDistribution;
}


BasicLog *EseConfiguration::GetReportLog() {
	if (logReport == NULL) 
		logReport = new AsciiLog(reportName.c_str(), extReport.c_str(), 0);
	return logReport;
}


BasicLog *EseConfiguration::GetPedigreeLog() {
	if (logPedigree == NULL) 
		logPedigree = new AsciiLog(reportName.c_str(), extPedigree.c_str(), 0);
	return logPedigree;
}


bool EseConfiguration::Validate() {
	bool success = true;	
	
/*	if (analysisStyle != EseConfiguration::PDT) {
		errorMsg = "This application is set to work with PDT analysis for now. The other analyses are still under construction\n";
		success=false;
	}
	else*/ 	if (analysisStyle == EseConfiguration::PDT && inputFormat != EseConfiguration::Pedigree) { 
		errorMsg = "PDT Requires the input format to be Pedigree\n";
		success=false;
	}
	else if (crossValidationCount < 1) {
		errorMsg = "Execution requires at least one for cross validation\n";
		success=false;
	}
	else if (bestDifferenceThreshold==0 && analysisStyle==MaxDiff) {
		errorMsg="Unable to perform Maximum Differential Analysis for threshold of 0! Please provide the MD metric\n";
		success=false;
	}
	else if (inputFile == "") {
		errorMsg="No input file provided.\n";
		success=false;
	}
	return success;
}


void EseConfiguration::ResetDistributions() {
	if (fDist)
		delete fDist;
	fDist = NULL;

	if (oDist)
		delete oDist;
	oDist = NULL;

	if (eDist)
		delete eDist;
	eDist = NULL;
}

/**
 * This MUST be called after the input has been parsed
 */
SnpEvalSuite *EseConfiguration::GetPdtEval(bool doRandomizeStatus, uint testNumber) {
	CaseControlStatus stat;
	BasicLog *report = NULL;
	//We don't want chatter during the P-Test runs
	if (!doRandomizeStatus)
		report = GetReportLog();

	EvalBalancedAccuracyPDT *evalPDT=new EvalBalancedAccuracyPDT(comboStart, comboEnd, minDThresh, doRandomizeStatus, report, doCalcThreshold, testNumber);

	//The status associated with this is a little different from the others. We need a vector
	GtLineParserMdrPdt *p=(GtLineParserMdrPdt *)((GtFileParserBuffered *)inputParser)->GetLineParser();


	FoldType pgStatus = NULL;


	if (doRandomizeStatus)
		p->SetupChildren(doRandomizeStatus);
	pgStatus = p->GetPgStatArray();	
		
	PdtFold *folds=p->GetFolds();
	p->GetStatusMask(&stat);
	evalPDT->SetFolds(folds, crossValidationCount, pgStatus, p->GetPgStatArrayCount());
	delete[] folds;
	//Setup the overall status
	evalPDT->SetOverallStatus(stat);
	//Setup the trainer (for x-fold validation later on)

	if (pTestCount > 0 && !lrDist)
		lrDist = new Omnibus(pTestCount, pTestSeed, "LR");
	if (pTestCount > 0 && !fDist)
		fDist = new NTest(comboEnd, pTestSeed, "MDR-PDT T");
	if (pTestCount > 0 && !oDist)
		oDist = new Omnibus(pTestCount, pTestSeed, "M. Odds Ratio");
	if (pTestCount > 0 && !eDist)
		eDist = new InvertedOmnibus(pTestCount, pTestSeed, "Pred. Error", 101.0);
		
	SnpEvalSuite *evalSuite=new SnpEvalSuite(evalPDT, fDist, oDist, eDist, lrDist);
		
	
	return evalSuite;
}

SnpEvalSuite *EseConfiguration::GetEval(CaseControlStatus &stat, GtFileParser *fileParser) {
	if (evalSuite)
		return evalSuite;

	evalSuite=new SnpEvalSuite();

	switch (analysisStyle) {
		case EseConfiguration::NoAnalysis:
			break;
		case EseConfiguration::MaxDiff:
		{
//			cout<<"Maximum Difference is still under consideration\n";
//			abort();
			//BalancedAccuracy maxdiff;
			EvalMaxDifference *md = new MDR::EvalMaxDifference(bestDifferenceThreshold, comboEnd);

			md->AppendStatus(stat);
			//md->status=stat;
			evalSuite->SetTrainer(md);
			
			//Build PTests
			fDist=new Omnibus(pTestCount, pTestCount, "Fitness Dist");
			evalSuite->SetFitnessDist(fDist);
/*			
			for (uint i=0; i<pTestCount; i++) {
				CaseControlStatus newStatus; 
				fileParser->AppendScrambledStatus(&newStatus);
				md = new ESE::EvalMaxDifference(bestDifferenceThreshold, comboEnd);
				md->AppendStatus(newStatus);
				evalSuite->AddPTest(md);
			}*/
			break;
		}
		case EseConfiguration::BalancedAccuracy:
		{
//			cout<<"Balanced Accuracy is still under consideration\n";
//			abort();
			EvalBalancedAccuracy *ba=new EvalBalancedAccuracy(comboEnd);

			//In the future, x-validation statuses will be added here. That will
			//Be the completed masks
			ba->AppendStatus(stat);
			//ba->status = stat;
			evalSuite->SetTrainer(ba);

			//Build PTests
			fDist=new Omnibus(pTestCount, pTestCount, "Fitness Dist");
			evalSuite->SetFitnessDist(fDist);
/*
			for (uint i=0; i<pTestCount; i++) {
				StatusContainer newStatus; 
				fileParser->AppendScrambledStatus(&newStatus);
				ba=new EvalBalancedAccuracy(comboEnd);
				total=newStatus.CombinedStatus();
				ba->AppendStatus(total);
				//eval->testSuite.push_back(newStatus);

				evalSuite->AddPTest(ba);
			}*/
			break;
		}
		case EseConfiguration::PDT:
		{
			cout<<"Internal Error: Inappropriate Stat for GetEval Call\n";
			abort();
			break;
		}


	}
	return evalSuite;
}

void EseConfiguration::SetFilename(const char *filename) {
	inputFile = filename;
	if (inputParser)
		delete inputParser;
	inputParser = NULL;

	if (families) {
		families->Purge();
		families->Release();
	}

}

void EseConfiguration::ValidateSnpCount( uint maxSnps) {
	if (comboStart > maxSnps) {
		comboStart = maxSnps;
		cout<<"** COMBO_START was set too high ("<<comboStart<<") for number of SNPs found in data. Set to "<<maxSnps<<" to fit the real data\n\n";
	}
	if (comboEnd > maxSnps) {
		comboEnd = maxSnps;
		cout<<"** COMBO_END was set too high ("<<comboEnd<<") for the number of SNPs found in data. Set to "<<maxSnps<<" to fit the real data\n\n";
	}
}

GtFileParser *EseConfiguration::GetInputFileParser(bool excludeLoci) {

	if (inputParser)
		return inputParser;
	
	if (inputFormat == MdrFormat) {
		GtLineParserMDR *lineParser=new GtLineParserMDR(1, false, maxGenotypes + 1 );
		inputParser=new GtFileParserBuffered(lineParser, inputFile.c_str());
	}
	else if (inputFormat == Pedigree) {
		families = 	FamilyRepository::Instance();
		BasicLog *pedLog = GetPedigreeLog();
		families->SetPedLog( pedLog );

		GtLineParserMdrPdt *lineParser=new GtLineParserMdrPdt(2, false, families, crossValidationCount, maxGenotypes +1 );

		if (pedExclusions.size() > 0)
			lineParser->ExcludePedigrees(pedExclusions);
		if (excludeLoci && lociExclusions.size() > 0)
			lineParser->ExcludeLoci(lociExclusions);

		lineParser->SetPedigreeStream(pedLog->GetStream());
		inputParser=new GtFileParserBuffered(lineParser, inputFile.c_str());

		mendelianErrorSearch = FixGenotypeErrors(mendelianSettings.actionLevel, mendelianSettings.pedThreshold, pedLog->GetStream());

		if (mendelianSettings.actionLevel > 0) {
			lineParser->AddPostEvaluation( &mendelianErrorSearch );
			
			if (mendelianSettings.actionLevel == 3) {
				mendelianStripper = StripBadLoci(&mendelianErrorSearch.affectedCounts, mendelianSettings.actionLevel, mendelianSettings.locThreshold);
				lineParser->AddPostEvaluation( &mendelianStripper );
			}
	
			if (doWriteCleanedDatafile && GetCleanedDataFile() )  {
				familyWriter = WriteFamilies(GetCleanedDataFile(), false);
				lineParser->AddPostEvaluation( &familyWriter );
			}
		}
	}

	else {	
		cerr<<"No handler for EseConfiguration::GetInputFileParser inputFormat=="<<inputFormat<<"\n";
		assert(0);
	}
	return inputParser;
}

void EseConfiguration::GenerateReport(ostream& os) {

	os<<"Configuration:\n";
	if (inputParser)
		inputParser->GenerateReport(os);

	os<<setw(45)<<right<<"Results Log: "<<GetReportLog()->GetFilename()<<endl;
	os<<setw(45)<<right<<"Distribution Log: "<<GetDistributionLog()->GetFilename()<<endl;
	if (analysisStyle) 
		os<<setw(45)<<right<<"Pedigree Log: "<<GetPedigreeLog()->GetFilename()<<endl;

	os<<setw(45)<<right<<"Model Size: "<<comboStart<<"-"<<comboEnd<<endl;

	if (crossValidationCount > 1)
    	os<<setw(45)<<"Cross Validation Folds: "<<crossValidationCount<<"\n";

	os<<setw(45)<<"Permutation Tests: "<<pTestCount<<"\n";
	os<<setw(45)<<"Permutation Test Seed: "<<pTestSeed<<"\n";
	os<<setw(45)<<right<<"Analysis Style: ";
	switch (analysisStyle) {
		case NoAnalysis:	
			os<<"No Analysis\n";
			break;
		case EseConfiguration::MaxDiff:
			os<<setw(45)<<right<<"Risk Assesment: ";
			if (doCalcThreshold)
				os<<"Based on each model's local threshold\n";
			else
				os<<"Based on overall Affected / Unaffected Ratio\n";
			cout<<"Maximum Difference Metric Single Snp. MD = ("<<bestDifferenceThreshold<<")\n";
			cout<<"\tOnly "<<comboEnd<<" SNP models will be reported on\n";
			if (reportModelCount >0)
				os<<"   * "<<reportModelCount<<" models will be reported (sorted by metric)\n";
			else {
				os<<"   * All models will be reported " ;
				if (doSortReport)
					os<<"(sorted by metric)\n";
				else
					os<<" unsorted.\n";
			}
			break;
		case EseConfiguration::BalancedAccuracy:
			os<<setw(45)<<right<<"Risk Assesment: ";
			if (doCalcThreshold)
				os<<"Based on each model's local threshold\n";
			else
				os<<"Based on overall Affected / Unaffected Ratio\n";
			break;
		case EseConfiguration::PDT:
			cout<<"Pedigree Analysis using PDT statistic\n";
			cout<<setw(45)<<right<<"Risk Assessment: "<<"Based on 1.0 (there is always a matched set of affecteds/unaffecteds)\n";
#ifdef USE_DOPTIMIZATION
			cout<<setw(45)<<" "<<"* D Optimization Set to: "<<minDThresh<<" (0.00 - 0.90)\n";
			cout<<setw(45)<<" "<<"  Higher values increase speed significantly but risk losing some models\n";
#endif
			break; 
		default:
			cout<<"Undefined Analysis Style ("<<analysisStyle<<")\n";
	}	
	os<<setw(45)<<"------------------------------------------"<<endl;
	os<<setw(40)<<"Mendelian Errors: "<<endl;

	if (mendelianSettings.actionLevel>0) {
		mendelianErrorSearch.Report(&os);
		mendelianStripper.Report(&os);
	}
	
	if (logLocus)
		os<<logLocus->ReportConfiguration();

}

}
