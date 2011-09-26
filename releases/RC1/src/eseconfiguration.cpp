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
#include "evalbalancedaccuracypdt.h"
#include "eseconfiguration.h"
#include "genetics/omnibusdistribution.h"
#include "genetics/ntestdistribution.h"
#include "matchedoddsratio.h"

//using namespace Genetics::Parser;
using namespace Genetics::Reporting;

namespace MDR {


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
char *EseConfiguration::ConfigFileValues::PedExclusionList			="PEDIGREE_EXCLUSIONS";
char *EseConfiguration::ConfigFileValues::VerboseMOR				="VERBOSE_MATCHED_ODDS_RATIO";
char *EseConfiguration::ConfigFileValues::VerbosePDT				="VERBOSE_EVALUATION";
char *EseConfiguration::ConfigFileValues::WriteFoldDetails			= "VERBOSE_FOLDING";


uint EseConfiguration::_instanceCount 								= 0;
EseConfiguration *EseConfiguration::_instance						= NULL;

EseConfiguration::EseConfiguration() : pTestCount(1000), pValueThreshold(0.30), minDThresh(0.0), pTestSeed(13), logLocus(NULL), comboEnd(1),comboStart(1), distribution(NULL),  evalSuite(NULL), inputParser(NULL), logDistribution(NULL), logReport(NULL), logPedigree(NULL), performPerfectSearch(false)
{
	families = 	NULL;
	errorMsg="";
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
	if (distribution)
		delete distribution;

	//Delete the logs
	if (logDistribution)
		delete logDistribution;
	if (logReport)
		delete logReport;
	if (logPedigree) 
		delete logPedigree;
	if (families)
		FamilyRepository::Release();
}

void EseConfiguration::SetDefaultValues() {
	reportName 				= "";
	extDistribution 		= "pdist";
	extLocus 				= "locus";
	extPedigree 			= "pedigree";
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
	else if (strcmp(key, ConfigFileValues::MinD)==0)
		minDThresh=atof(value);
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
		doCreateVirtualSibs=GetBoolean(value);
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
		LoadList(pedExclusions, remainder);
	else if (strcmp(key, ConfigFileValues::LociExclusionList) == 0)
		LoadList(lociExclusions, remainder);
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
	else
		success=false;
	return success;
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
	
	if (analysisStyle != EseConfiguration::PDT) {
		errorMsg = "This application is set to work with PDT analysis for now. The other analyses are still under construction\n";
		success=false;
	}
	else 	if (analysisStyle == EseConfiguration::PDT && inputFormat != EseConfiguration::Pedigree) { 
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

/**
 * This MUST be called after the input has been parsed
 */
SnpEvalSuite *EseConfiguration::GetPdtEval(bool doRandomizeStatus, uint testNumber) {
	CaseControlStatus stat;
	BasicLog *report = NULL;
	//We don't want chatter during the P-Test runs
	if (!doRandomizeStatus)
		report = GetReportLog();

	SnpEvalSuite *evalSuite=new SnpEvalSuite();
	EvalBalancedAccuracyPDT *evalPDT=new EvalBalancedAccuracyPDT(comboStart, comboEnd, minDThresh, doRandomizeStatus, report, doCalcThreshold, testNumber);

	//The status associated with this is a little different from the others. We need a vector
	GtLineParserMdrPdt *p=(GtLineParserMdrPdt *)((GtFileParserBuffered *)inputParser)->GetLineParser();

	FoldType pgStatus = NULL;

	if (doRandomizeStatus)
		p->SetupChildren(doRandomizeStatus);
	else 
		pgStatus = p->GetPgStatArray();	
		
	PdtFold *folds=p->GetFolds();
	p->GetStatusMask(&stat);
	evalPDT->SetFolds(folds, crossValidationCount, pgStatus, p->GetPgStatArrayCount());
	
	//Setup the overall status
	evalPDT->SetOverallStatus(stat);
	//Setup the trainer (for x-fold validation later on)
	evalSuite->SetTrainer(evalPDT);

	//Build PTests
	if (distribution == NULL)
		distribution=new NTestDistribution(pTestSeed, comboEnd, pTestCount, GetDistributionLog());
	evalSuite->SetDistribution(distribution);
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
			break;
		}
		case EseConfiguration::BalancedAccuracy:
		{
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



GtFileParser *EseConfiguration::GetInputFileParser(bool doExcludeLoci) {

	if (inputParser)
		return inputParser;
	
	if (inputFormat == Pedigree) {
		families = 	FamilyRepository::Instance();
		BasicLog *pedLog = GetPedigreeLog();
		families->SetPedLog( pedLog );

		GtLineParserMdrPdt *lineParser=new GtLineParserMdrPdt(2, false, families, crossValidationCount, maxGenotypes +1 );
		if (pedExclusions.size() > 0)
			lineParser->ExcludePedigrees(pedExclusions);
		if (doExcludeLoci && lociExclusions.size() > 0)
			lineParser->ExcludeLoci(lociExclusions);

		lineParser->DoCreateVirtualSibs(doCreateVirtualSibs);
		lineParser->SetPedigreeStream(pedLog->GetStream());
		inputParser=new GtFileParserBuffered(lineParser, inputFile.c_str());
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
	else	{
		os<<"Configuration Error: Unable to acquire input parser.\n";
	}

	os<<setw(45)<<right<<"Results Log: "<<GetReportLog()->GetFilename()<<endl;
	os<<setw(45)<<right<<"Distribution Log: "<<GetDistributionLog()->GetFilename()<<endl;
	if (analysisStyle) 
		os<<setw(45)<<right<<"Pedigree Log: "<<GetPedigreeLog()->GetFilename()<<endl;

	os<<setw(45)<<right<<"Model Size: "<<comboStart<<"-"<<comboEnd<<endl;




	os<<setw(45)<<"Permutation Tests: "<<pTestCount<<"\n";
	os<<setw(45)<<"Permutation Test Seed: "<<pTestSeed<<"\n";
	os<<setw(45)<<right<<"Analysis Style: ";
	switch (analysisStyle) {
		case NoAnalysis:	
			os<<"No Analysis\n";
			break;
		case EseConfiguration::PDT:
			cout<<"Pedigree Analysis using PDT statistic\n";
			cout<<setw(45)<<right<<"Risk Assessment: "<<"Based on 1.0 \n";
			break; 
		default:
			cout<<"Undefined Analysis Style ("<<analysisStyle<<")\n";
	}	
	if (logLocus)
		os<<logLocus->ReportConfiguration();

}

}
