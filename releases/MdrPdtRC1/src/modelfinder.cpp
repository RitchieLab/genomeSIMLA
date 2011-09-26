//
// C++ Implementation: eseapp
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "modelfinder.h"
#include "evalbalancedaccuracypdt.h"
#include <boost/timer.hpp>

#ifndef TEST_APP
int main(int argc, char **argv) {
	MDRPDT::ModelFinder app;

	if (app.ParseCmdLine(argc, argv))
		app.Start();
	return 0;
}
	
#endif

namespace MDRPDT {

using namespace boost;

ModelFinder::ModelFinder() : doAnalyze(true) {
	appname="mdr-pdt";
	appfunction="Multifactor Dimensionality Reduction - Pedigree Disequilibrium Test";
	authors="Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu) & Eric Torstenson (torstenson@chgr.mc.vanderbilt.edu)";
	major=1;
	minor=0;
	configuration = AppConfiguration::Instance();
}


ModelFinder::~ModelFinder()
{
	AppConfiguration::Release();
}


bool ModelFinder::DumpSampleConfig(int argc, char **argv) {
	bool handled = false;
	stringstream ss;
	stringstream specialExt;
	ss<<"/******************** Basic Settings\n";
	ss<<" * Describe the type of models of interest. Models consist on 1 or more SNPs.\n";
	ss<<" * The minimum number of SNPs in a model to be investigated. Valid settings: 1..COMBO_END\n";
	string comboStart = "1";
	string comboEnd = "2";
	if (argc > 3)
		comboStart = argv[3];
	if (argc > 4)
		comboEnd  = argv[4];
	ss<<"COMBO_START            "<<comboStart<<"\n";
	ss<<" *  The maximum number of SNPs to be considered. Valid settings: COMBO_START...\n";
	ss<<"COMBO_END              "<<comboEnd<<"\n";
	ss<<" *  If the data has little or no missing data, it will run faster by using the overall threshold\n";
	ss<<"USE_MODEL_THRESHOLD    ON\n";

	//Produce a PDT formatted configuration file
	if (strcmp(argv[1], AppConfiguration::ConfigFileValues::PDTAnalysis) == 0) {
		handled = true;
		ss<<" *  The value used to indicate that an individual is affected\n";
		ss<<"AFFECTED_VALUE         2\n";
		ss<<" *  The value used to indicate that an individual is unaffected\n";
		ss<<"UNAFFECTED_VALUE       1\n";
		ss<<" *  All other individuals will be considered to be of unknown status and will not contribute to the calculations\n";
		ss<<"\n\n/****************** Input format\n";
		ss<<" * PEDIGREE   - For PDT anayses there is only one format supported\n";
		ss<<"INPUTFORMAT            PEDIGREE\n";
		ss<<" *  The name of the input file where your SNP data is to be found. This file must be space delimited\n";
		string filename = "YourFilename.ped";
		if (argc>2)
			filename=argv[2];
		ss<<"INPUTFILE              "<<filename<<"\n";
		specialExt<<"EXT_PEDIGREE          pedigree      //Pedigree related errors and notes are logged here\n";
		ss<<" *  There is only one Pedigree analysis currently available\n";
		ss<<"ANALYSISSTYLE          PDT\n";
		ss<<" *  Turn On/Off Triad expansion\n";
		ss<<"EXPAND_TRIOS           Yes\n";
#ifdef USE_DOPTIMIZATION
		ss<<" *  Experimental Optimization. Valid Optimization settings (0.1 - 0.9). \n";
		ss<<"DOPTIMIZATION			0.10\n";
#endif 
		ss<<" *  Matched Odds Verbose On/Off. Turning this on generates the table of pair contributions\n";
		ss<<"VERBOSE_MATCHED_ODDS_RATIO Off\n";
		ss<<" *  Evaluation Verbose Yes/No\n";
		ss<<"VERBOSE_EVALUATION     No\n";
	}



	ss<<"\n\n/****************** Reporting\n";
	ss<<" * Reports are set up in the form: REPORT_NAME.EXT Each type or report has it's own different extension.\n";
	ss<<" * Any report whose extension is STDOUT will be redirected to standard out instead of written to a file.\n";
	ss<<" * Any report whose extension is NOLOG will be skipped altogether.\n";
	ss<<" * By default, REPORT_NAME is just the name of the configuration file\n";  
	ss<<" *REPORT_NAME MyReportName\n";  
	ss<<"EXT_DISTRIBUTION       pdist      //This is where the distribution of the p-tests is written\n";
	ss<<"EXT_REPORT             STDOUT     //This is where the final results are written\n";
	ss<<specialExt.str();

	ss<<"\n\n/****************** Permutation Tests\n";
	ss<<" * Number of permutation runs to be executed. 1000 is recommended.\n";
	ss<<"PTEST_COUNT            1000\n";
	ss<<" * The seed associated with the tests. Each test gets a new seed\n";
	ss<<"PTEST_SEED             333\n";
	ss<<" * Reporting can be done based on a p-value threshold. Any model whose significance exceeds this value will be reported.\n";
	ss<<"PVAL_THRESHOLD         0.05\n";
	
	if (handled)
		cout<<ss.str();
	return handled;
}

bool ModelFinder::ParseCmdLine(int argc, char **argv) {
	bool success=false;
	if (argc < 2) 
		PrintHelp();
	else {
		if (!DumpSampleConfig(argc, argv)) {
			configFile = argv[1];
	
			LineParser lp;
			lp.Parse(configFile.c_str(), configuration);
			configuration->reportName = configFile;
			success=configuration->Validate();
			
			if (!success)
				cout<<configuration->errorMsg;
			else 
				configuration->PostLoad();
			
			for (int i =2; i < argc; i++) {
				doAnalyze=false;
				models.push_back(argv[i]);
			}
		}
	}
	return success;
}


bool ModelFinder::VerifyConfiguration() {
	bool isGood = true;

	return isGood;
}


/**
 * Prints the help contents
 */
void ModelFinder::PrintHelp()	{
	PrintBanner();
	cout<<"usage: "<<appname<<" <configuration file> [optional model IDs]\n";
	cout<<"\nThis search tool can be used to perform analyses on SNP data or to explore a single model\n";
	cout<<"* To execute analysis, simply run the command with the configuration file only\n";
	cout<<"* To explore a specific model, specify the configuration file and one or more models separated by space\n";
	cout<<"* To dump an example configuration file, pass the name of the analyses desired as the sole paramter. Valid choices are [PDT]\n"; 
}

/**	
 * Starts execution
 */
void ModelFinder::Start()	{
	PrintBanner();

	//Pedigree searches behave very differently- gotta switch approach accordingly
	if (configuration->analysisStyle == AppConfiguration::PDT)
		if (doAnalyze)
			PedigreeSearch();
		else
			PedigreeDescribeModels();
	else
		cout<<"Invalid Search method selected. Currently, options are restricted to PDT searches on pedigree data.\n";
}

void ModelFinder::PedigreeDescribeModels() {
	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser();
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();

	snps.ParseInputFile( fileParser, NULL );

	//Perform some basic reporting
	configuration->GenerateReport(cout);
	
	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	EvalBalancedAccuracyPDT *evalPDT=new EvalBalancedAccuracyPDT( 1, 1, 0.0, false, NULL, configuration->doCalcThreshold);
	FoldType st=pdtParser->GetStatArray();

	FoldType pgStatus = pdtParser->GetPgStatArray();

	evalPDT->AppendStatus(st, pdtParser->GetStatArrayCount(), pgStatus, pdtParser->GetPgStatArrayCount());
	evalPDT->SetOverallStatus(stat);

	//Get the number of models we are about to explore
	uint count = models.size();
	SnpAligned *model;
	string modelID;
	for (uint i=0; i<count; i++) {
		modelID=models[i];
		model = snps.GetSnp(modelID.c_str());
		if (model)	{
			//Evaluate the model with with high output
			evalPDT->EvaluateVerbose(model);
			if (i+1<count)
				cout<<"\n\n-------------------------------------------------------------------------------\n\n";
			model->ReduceInstanceCount();
		}

	}
}
void ModelFinder::PedigreeSearch() {

	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser();

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );

	if (snps.GetSnpCount() == 0) {
		cout<<"\nNo data to search over. Please check that the file that you specified exists\n";
		abort();
	}
	timer progress;

	//Perform some basic reporting
	configuration->GenerateReport(cout);

	//Close the locus log, if it exists
	if (configuration->logLocus)
		configuration->logLocus->Close();

	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	//Create a report to store the top models
	SnpReposTxtSorted bestScore(configuration->reportModelCount);
	SnpEvalSuite *e = configuration->GetPdtEval( false, 0 );

	//The distribution object needs to be passed to the various ptests
	PermutationTestDist *dist = e->GetDistribution();
	bestScore.distribution = dist;
	BasicLog *reportLog = configuration->GetReportLog();

	ostream *sigReport = reportLog->GetStream();

	//If some analysis is required, do it now
	if (e) {
		progress.restart();
		//Perform the search
		if (EvalBalancedAccuracyPDT::Verbose)
			cout<<setw(16)<<right<<"Model"<<setw(8)<<"   T Statistic"<<endl;
		else
			cout<<"  Searching\n";

		uint modelCount = snps.Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, e);
		cout<<"\n\nPedigree Search of "<<e->GetSnpCount()<<" : "<<modelCount<<" completed in "<<progress.elapsed()<<" seconds\n";

		
		cout<<"\n\nPerforming "<<configuration->pTestCount<<" Permutation Tests:\n";
		
		//GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();
		cout<<setw(6)<<"Test #"<<setw(15)<<"Load Time (s)"<<setw(30)<<"Execution Time(s)"<<setw(15)<<"Models Seen"<<endl;
		for (uint i=0; i<configuration->pTestCount; i++) {
			progress.restart();
			SnpRepository ptest;

			//This is a quick fix for random numbers. 
			srand(configuration->pTestSeed + i);

			//Acquire the evaluation suite- this time, we want randomized status
			SnpEvalSuite *pEval = configuration->GetPdtEval(true, i);

			//Populate the test repository
			ptest.LoadData( fileParser);
			//Acquire the status. Technically, this should not change, but, in case the order changed for some reason, we want to make 
			//sure we are up to date.
			fileParser->GetStatusMask(&stat);

			cout<<setw(6)<<i<<setw(15)<<setprecision(3)<<progress.elapsed();
			progress.restart();
			//Perform search on test data
			uint modelCount = ptest.Evaluate( configuration->comboStart - 1, configuration->comboEnd - 1, NULL, pEval);
			cout<<setw(30)<<progress.elapsed()<<setw(15)<<modelCount<<"\n";
			delete pEval;
		}
		//Set up the distribution and sort the results
		e->BuildDistribution();

		//Report on the distribution
		if (e->GetDistribution())
			e->GetDistribution()->Report();
	
		float minPVal = configuration->pValueThreshold;	
		
		bestScore.Sort();
		if (sigReport == NULL) {
			cout<<"\n\n";
			e->ReportResults(&cout);
		}
		else {
			//Eventually, this will go into a specialized report
			for (uint modelSize = configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
				uint padding = 4 * modelSize + 4;
				uint  reportSize=bestScore.GetTotalEntryCount(modelSize);
				*sigReport<<"\n\nDistribution Report "<<modelSize + 1<<" SNP(s) per model"<<endl;
				*sigReport<<"\t* Only models with a p-value < "<<setprecision(3)<<minPVal<<" are reported. "<<endl;
				*sigReport<<setw(padding)<<""<<setw(12)<<"T    "<<setw(16)<<"Statistical"<<endl;
				*sigReport<<setw(padding)<<"Model"<<setw(12)<<"Statistic"<<setw(16)<<"Significance"<<endl;
				
				float pValue=1.0;
				SnpReposTxtSorted::ReportModel *model = bestScore.GetReportModel(modelSize, 0);
				if (model) 
					pValue=dist->GetPValue(modelSize+1, model->score);
	
				bool significantFinds=false;
				for (uint idx=1; idx < reportSize; idx++) {
					if (pValue <= minPVal) {
						significantFinds=true;
						*sigReport<<setw(padding)<<model->GetLabel();
						*sigReport<<setw(12)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3)<<model->score;
						*sigReport<<"       p < "<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3)<<pValue<<endl;
					}
					model = bestScore.GetReportModel(modelSize, idx);
					pValue=dist->GetPValue(modelSize+1, model->score);
					
				}
				if (!significantFinds)
					*sigReport<<setw(12)<<"None found"<<endl;
			}
	
			//Report on the overal resutlts
			*sigReport<<"\n\n";
			e->ReportResults(sigReport);
		}
		delete e;

	}
}




}
