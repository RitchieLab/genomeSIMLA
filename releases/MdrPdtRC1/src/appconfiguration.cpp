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
#include "appconfiguration.h"
#include "genetics/omnibusdistribution.h"
#include "genetics/ntestdistribution.h"
#include "genetics/gtfileparserbuffered.h"
#include "matchedoddsratio.h"
#include "gtlineparsermdrpdt.h"

//using namespace Genetics::Parser;
using namespace Genetics::Reporting;

namespace MDRPDT {


////Fundamental Settings
									//The first locus count for models of interest
char *AppConfiguration::ConfigFileValues::ComboStart				= "COMBO_START";
									//The largest size of interest
char *AppConfiguration::ConfigFileValues::ComboEnd					= "COMBO_END";			
char *AppConfiguration::ConfigFileValues::CalcThreshold				= "USE_MODEL_THRESHOLD";
char *AppConfiguration::ConfigFileValues::AffectedValue 			= "AFFECTED_VALUE";
char *AppConfiguration::ConfigFileValues::UnaffectedValue			= "UNAFFECTED_VALUE";

////Settings for overrideing various report extensions
									//Let user set the extension for distribution reports
char *AppConfiguration::ConfigFileValues::ExtDistribution			= "EXT_DISTRIBUTION";
									//Let user set the extension for the overall report
char *AppConfiguration::ConfigFileValues::ExtReport					= "EXT_REPORT";
									//Let user set the extension used for locus reports		
char *AppConfiguration::ConfigFileValues::ExtLocus					= "EXT_LOCUS";
									//Let user set the extension used for the pedigree report
char *AppConfiguration::ConfigFileValues::ExtPedigree				= "EXT_PEDIGREE";
									//Let user set the report name
char *AppConfiguration::ConfigFileValues::ReportName				= "REPORT_NAME";					

////Overriding other reporting details
								//Cutoff for the PValue output
char *AppConfiguration::ConfigFileValues::MinPVal					= "PVAL_THRESHOLD";		
								//Number of PTests
char *AppConfiguration::ConfigFileValues::PTestCount				="PTEST_COUNT";
								//The seed used for the runs
char *AppConfiguration::ConfigFileValues::PTestSeed					="PTEST_SEED";
	
		//Overriding input details
									//<Used to override the default input format
char *AppConfiguration::ConfigFileValues::_InputFile				="INPUTFILE";
									//Used to override the default format of input data
char *AppConfiguration::ConfigFileValues::_InputFormat				="INPUTFORMAT";
									//<Used to specify to use MDR formatted input
char *AppConfiguration::ConfigFileValues::_MdrFormat				="MDRFORMAT";
									//Used to specify pedigree formatted input
char *AppConfiguration::ConfigFileValues::PedigreeFormat			="PEDIGREE";

		//Overriding how the data is to be interpretted
									//<Use to override the max genotypes
char *AppConfiguration::ConfigFileValues::_GenotypeCount			="GENOTYPES";
									//<Used to override the type of analysis
char *AppConfiguration::ConfigFileValues::_AnalysisStyle			="ANALYSISSTYLE";
									//<Used to specify Max Difference on CC data			
char *AppConfiguration::ConfigFileValues::_MaxDiff					="MAXIMIZEDIFF";
									//<Used to specify Balanced Accuracy Analyses on CC data
char *AppConfiguration::ConfigFileValues::BalancedAccuracy 			="BALANCEDACCURACY";
									//<Used to specify PDT analyses on pedigree data
char *AppConfiguration::ConfigFileValues::PDTAnalysis				="PDT";
		
//MD Specific Settings
									//Set the number of models reported on
char *AppConfiguration::ConfigFileValues::_ReportModelCount			="REPORTMODELCOUNT\0";
									//Set the number of affected (required for Base Pair layout)
char *AppConfiguration::ConfigFileValues::_AffectedCount			="AFFECTEDCOUNT\0";
									//Set the number of unaffected (required for Base Pair layout)
char *AppConfiguration::ConfigFileValues::_UnaffectedCount			="UNAFFECTEDCOUNT\0";
									//Set the difference threshold (minimum required to be reported on)
char *AppConfiguration::ConfigFileValues::_BestDifferenceThreshold	="BESTDIFFERENCETHRESHOLD\0";
									//set whether or not the report is sorted
char *AppConfiguration::ConfigFileValues::_SortReportModels			="PERFORMSORT\0";
	
//Pedigree Specific settings
									//<The percentage within the highest D to be considered
char *AppConfiguration::ConfigFileValues::MinD						="DOPTIMIZATION";
char *AppConfiguration::ConfigFileValues::CreateVirtualSibs			="EXPAND_TRIOS";
char *AppConfiguration::ConfigFileValues::VerboseMOR				="VERBOSE_MATCHED_ODDS_RATIO";
char *AppConfiguration::ConfigFileValues::VerbosePDT				="VERBOSE_EVALUATION";


uint AppConfiguration::_instanceCount 								= 0;
AppConfiguration *AppConfiguration::_instance						= NULL;

AppConfiguration::AppConfiguration() : pTestCount(1000), pValueThreshold(0.30), minDThresh(0.0), pTestSeed(13), logLocus(NULL), comboEnd(1),comboStart(1), distribution(NULL),  evalSuite(NULL), inputParser(NULL), logDistribution(NULL), logReport(NULL), logPedigree(NULL), performPerfectSearch(false)
{
	families = 	NULL;
	errorMsg="";
	SetDefaultValues();
	PostLoad();
}


AppConfiguration::~AppConfiguration()
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

void AppConfiguration::SetDefaultValues() {
	reportName 				= "";
	extDistribution 		= "pdist";
	extLocus 				= "locus";
	extPedigree 			= "pedigree";
	extReport 				= "STDOUT";

	locusWriteStyle			= LocusLogAscii::WriteAll;
	analysisStyle			= MaxDiff;
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
}

void AppConfiguration::PostLoad() {
	if (analysisStyle != PDT && logLocus) {
		delete logLocus;
		logLocus = new LocusLogAscii(locusWriteStyle, reportName.c_str(), extPedigree.c_str());
	}
}

void AppConfiguration::SetInputFormat(const char *value) {
	if (strcmp(value, ConfigFileValues::_MdrFormat)==0) 
		inputFormat=MdrFormat;
	else if (strcmp(value, ConfigFileValues::PedigreeFormat)==0)
		inputFormat=Pedigree;

}

void AppConfiguration::SetAnalysisStyle(const char *value) {
	if (strcmp(value, ConfigFileValues::BalancedAccuracy)==0)
		analysisStyle=BalancedAccuracy;
	else if (strcmp(value, ConfigFileValues::PDTAnalysis)==0)
		analysisStyle=PDT;
	else if (strcmp(value, ConfigFileValues::_MaxDiff)==0) {
		analysisStyle=MaxDiff;
		performPerfectSearch = true;
	}

}



bool AppConfiguration::SetValues(const char *key, const char *value) {
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


BasicLog *AppConfiguration::GetDistributionLog() {
	if (logDistribution == NULL) 
		logDistribution = new AsciiLog(reportName.c_str(), extDistribution.c_str(), 0);
	return logDistribution;
}


BasicLog *AppConfiguration::GetReportLog() {
	if (logReport == NULL) 
		logReport = new AsciiLog(reportName.c_str(), extReport.c_str(), 0);
	return logReport;
}


BasicLog *AppConfiguration::GetPedigreeLog() {
	if (logPedigree == NULL) 
		logPedigree = new AsciiLog(reportName.c_str(), extPedigree.c_str(), 0);
	return logPedigree;
}


bool AppConfiguration::Validate() {
	bool success = true;	
	
/*	if (analysisStyle != AppConfiguration::PDT) {
		errorMsg = "This application is set to work with PDT analysis for now. The other analyses are still under construction\n";
		success=false;
	}
	else*/ if (analysisStyle == AppConfiguration::PDT && inputFormat != AppConfiguration::Pedigree) { 
		errorMsg = "PDT Requires the input format to be Pedigree\n";
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
SnpEvalSuite *AppConfiguration::GetPdtEval(bool doRandomizeStatus, uint testNumber) {
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
		
	FoldType st=p->GetStatArray();
	p->GetStatusMask(&stat);
	
	//Append the status array to the Evaluation object
	evalPDT->AppendStatus(st, p->GetStatArrayCount(), pgStatus, p->GetPgStatArrayCount());
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


GtFileParser *AppConfiguration::GetInputFileParser() {

	if (inputParser)
		return inputParser;
	else if (inputFormat == Pedigree) {
		families = 	FamilyRepository::Instance();
		BasicLog *pedLog = GetPedigreeLog();
		families->SetPedLog( pedLog );

		GtLineParserMdrPdt *lineParser=new GtLineParserMdrPdt(2, false, families, maxGenotypes +1 );
		lineParser->DoCreateVirtualSibs(doCreateVirtualSibs);
		lineParser->SetPedigreeStream(pedLog->GetStream());
		inputParser=new GtFileParserBuffered(lineParser, inputFile.c_str());
	}

	else {	
		cerr<<"No handler for AppConfiguration::GetInputFileParser inputFormat=="<<inputFormat<<"\n";
		assert(0);
	}
	return inputParser;
}

void AppConfiguration::GenerateReport(ostream& os) {

	os<<"Configuration:\n";
	if (inputParser)
		GetInputFileParser()->GenerateReport(os);
	else	{
		os<<"Configuration Error: Unable to acquire input parser.\n";
	}

	os<<setw(45)<<right<<"Results Log: "<<GetReportLog()->GetFilename()<<endl;
	os<<setw(45)<<right<<"Distribution Log: "<<GetDistributionLog()->GetFilename()<<endl;
	if (analysisStyle) 
		os<<setw(45)<<right<<"Pedigree Log: "<<GetPedigreeLog()->GetFilename()<<endl;

	os<<setw(45)<<right<<"Model Size: "<<comboStart<<"-"<<comboEnd<<endl;
	os<<setw(45)<<right<<"Risk Assesment: ";
	if (doCalcThreshold)
		os<<"Based on each model's local threshold\n";
	else
		os<<"Based on overall Affected / Unaffected Ratio\n";
	os<<setw(45)<<right<<"Analysis Style: ";
	switch (analysisStyle) {
		case NoAnalysis:	
			os<<"No Analysis\n";
			break;
		case AppConfiguration::PDT:
			cout<<"Pedigree Analysis using PDT statistic\n";
#ifdef USE_DOPTIMIZATION
			cout<<setw(45)<<" "<<"* D Optimization Set to: "<<minDThresh<<" (0.00 - 0.90)\n";
			cout<<setw(45)<<" "<<"  Higher values increase speed significantly but risk losing some models\n";
#endif
			break; 
		default:
			cout<<"Undefined Analysis Style ("<<analysisStyle<<")\n";
	}	
	if (logLocus)
		os<<logLocus->ReportConfiguration();

}

}
