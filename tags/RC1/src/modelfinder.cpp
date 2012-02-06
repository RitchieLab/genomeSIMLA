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


//This is required for the way the stdlib isn't supposed to use rand(). If we need to 
//generate an random integer, let's use these since we might have different random 
//generators running concurrently otherwise
#ifdef _GLIBCPP_HAVE_DRAND48
#define SRAND(A) srand48(A)
#define NEXTRAND() lrand48()
#else
#define SRAND(A) srand(A)
#define NEXTRAND() rand()
#endif


#ifndef TEST_APP
int main(int argc, char **argv) {
	MDR::ModelFinder app;

	if (app.ParseCmdLine(argc, argv))
		app.Start();
	return 0;
}
	
#endif

namespace MDR {

using namespace boost;

ModelFinder::ModelFinder() : doAnalyze(true) {
	appname="mdr-pdt";
	appfunction="Multifactor Dimensionality Reduction - Pedigree Disequilibrium Test";
	authors="Marylyn Ritchie & Eric Torstenson\nPlease forward any comments or errors to: mdr-pdt@chgr.mc.vanderbilt.edu";
	major=1;
	minor=1;
	bugFixes = 4;
	configuration = EseConfiguration::Instance();
}


ModelFinder::~ModelFinder()
{
	EseConfiguration::Release();
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




	ss<<" * You can exclude loci from analyses by adding them to the following variable\n";
	ss<<" * The application does not actually keep the data associated with these loci, but they \n";
	ss<<" * do retain the positions, so model 13x22 with none excluded would be the same as \n";
	ss<<" * 13x22 with 14, 15 and 16 excluded\n";
	ss<<"EXCLUDE_LOCUS	\n";
	ss<<" * The following would exclude loci 1, 4 and 31 from analyses\n";
	ss<<"*EXCLUDE_LOCUS 1 31 4\n";


	//Produce a PDT formatted configuration file
	if (strcmp(argv[1], EseConfiguration::ConfigFileValues::PDTAnalysis) == 0) {
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
		
		ss<<" *  You can exclude certain pedigrees from analysis. Buy removing the asterisk from the line below\n";
		ss<<" *  will tell the application to exclude pedigrees 13, 21 and 1003\n";
		ss<<"*PEDIGREE_EXCLUSIONS 13 21 1003\n";

		ss<<"ANALYSISSTYLE          PDT\n";
		ss<<" *  Turn On/Off Triad expansion\n";
		ss<<"EXPAND_TRIOS           Yes\n";

		ss<<" *  Evaluation Verbose Yes/No\n";
		ss<<"VERBOSE_EVALUATION     No\n";




	}

	if (strcmp(argv[1], EseConfiguration::ConfigFileValues::BalancedAccuracy) == 0) {
		handled = true;
		ss<<"\n\n/****************** Input Format\n";
		ss<<" * MDRFORMAT - For performing basic Case Control based SNP evaluation\n";
		ss<<"INPUTFORMAT             MDRFORMAT\n";
		ss<<" *  The name of the input file where your SNP data is to be found. This file must be space delimited\n";
		string filename = "YourFilename.mdr";
		if (argc>2)
			filename=argv[2];
		ss<<"INPUTFILE              "<<filename<<"\n";
		ss<<" *  There is only one format that for Case/Control analysis currently available\n";
		ss<<"ANALYSISSTYLE          "<<EseConfiguration::ConfigFileValues::BalancedAccuracy<<"\n";
	}

	if (strcmp(argv[1], EseConfiguration::ConfigFileValues::_MaxDiff) == 0) {
		handled = true;
		ss<<"\n\n/****************** Input Format\n";
		ss<<" * MDRFORMAT - For performing basic Case Control based SNP evaluation\n";
		ss<<"INPUTFORMAT             MDRFORMAT\n";
		ss<<" *  The name of the input file where your SNP data is to be found. This file must be space delimited\n";
		string filename = "YourFilename.ese";
		if (argc>2)
			filename=argv[2];
		ss<<"INPUTFILE              "<<filename<<"\n";
		ss<<" *  There is only one format that for Case/Control analysis currently available\n";
		ss<<"ANALYSISSTYLE          "<<EseConfiguration::ConfigFileValues::_MaxDiff<<"\n";
		ss<<" *  Allows you to set a limit on the number of models reported.\n";
		ss<<" *  When set to 0, all will be written to the screen\n";
		ss<<"REPORTMODELCOUNT        0\n";
		ss<<" *  Threshold at which point the application will consider a model to be of interest.\n";
		ss<<" *  This value will differ from data set to data set. If this is set too high, the application\n";
		ss<<" *  will still show the highest model seen and it's difference. \n";
		ss<<"BESTDIFFERENCETHRESHOLD 64\n";
		ss<<" *  Indicate if the application should sort the final report (only relevant if REPORTMODELCOUNT > 0\n";
		ss<<"PERFORMSORT             Yes\n";

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
	ss<<"PTEST_SEED             1371\n";
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
	if (configuration->analysisStyle == EseConfiguration::PDT)
		if (doAnalyze)
			PedigreeSearch();
		else
			PedigreeDescribeModels();
}

void ModelFinder::PedigreeDescribeModels() {
	SRAND(configuration->pTestSeed );

	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser(false);
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();

	snps.ParseInputFile( fileParser, NULL );

	//Perform some basic reporting
	configuration->GenerateReport(cout);
	
	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	EvalBalancedAccuracyPDT *evalPDT=new EvalBalancedAccuracyPDT( 1, 1, 0.0, false, NULL, configuration->doCalcThreshold);
	PdtFold *folds = pdtParser->GetFolds();
	FoldType pgStatus = pdtParser->GetPgStatArray();
	evalPDT->SetFolds(folds, configuration->crossValidationCount, pgStatus, pdtParser->GetPgStatArrayCount());
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
		else 
			cout<<"\n!! Model '"<<modelID<<"' returned an invalid model.Please see the User's Guide for information on formatting your model query properly. \n";

	}
}
void ModelFinder::PedigreeSearch() {
	SRAND(configuration->pTestSeed );

	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser(true);

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
			SRAND(configuration->pTestSeed + i + 1);

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
