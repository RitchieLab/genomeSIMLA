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
#ifdef USE_MPI
#include "mpi.h"
#include <pthread.h>
pthread_mutex_t input_lock  = PTHREAD_MUTEX_INITIALIZER;
#define TESTLOCK pthread_mutex_lock(&input_lock)
#define TESTUNLOCK pthread_mutex_unlock(&input_lock)

#else
//We don't want these to mean anything if we arne't using MPI
#define TESTLOCK
#define TESTUNLOCK

#endif 	//USE_MPI

#include "modelfinder.h"
#include "evalbalancedaccuracypdt.h"
//#include "snprepostxtsorted.h"
#include <boost/timer.hpp>




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

ModelFinder::ModelFinder() : doAnalyze(true), nodeID(0) {
	appname="mdr-pdt";
	appfunction="Multifactor Dimensionality Reduction - Pedigree Disequilibrium Test";
	authors="Marylyn Ritchie & Eric Torstenson\nPlease forward any comments or errors to: mdr-pdt@chgr.mc.vanderbilt.edu";
	major=1;
	minor=5;
	bugFixes=0;
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
#ifdef CROSS_VALIDATION
	ss<<" * Set the number of cross validation folds (we recommend 1, 5 or 10. 1 Is normal execution)\n";
	ss<<"CROSSVALINTERVAL       1\n";
#endif
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
		ss<<"EXPAND_TRIOS           On\n";
		/** We need to verify this first- need some data to do that with
		ss<<" *  Turn On/Off Expansion of all affected sibs (when parental data is present.)\n";
		ss<<"EXPAND_ALL             Off\n";
		*/
#ifdef USE_DOPTIMIZATION
		ss<<" *  Experimental Optimization. Valid Optimization settings (0.1 - 0.9). \n";
		ss<<"DOPTIMIZATION			0.10\n";
#endif 
		ss<<" *  Matched Odds Verbose On/Off. Turning this on generates the table of pair contributions\n";
		ss<<"VERBOSE_MATCHED_ODDS_RATIO Off\n";
		ss<<" *  Evaluation Verbose Yes/No\n";
		ss<<"VERBOSE_EVALUATION     No\n";
#ifdef CROSS_VALIDATION
		ss<<" *  Writes the arrangement of the pedigree assignments (for cross validation folding) to the pedigree report\n";
		ss<<"VERBOSE_FOLDING        No\n";
#endif 
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
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
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
	masterID = 0;

#ifdef USE_MPI
	nodeID = MPI::COMM_WORLD.Get_rank();
	processCount = MPI::COMM_WORLD.Get_size();
#endif

	if (nodeID == 0)
		PrintBanner();

	//Pedigree searches behave very differently- gotta switch approach accordingly
	if (configuration->analysisStyle == EseConfiguration::PDT)
		if (doAnalyze) {
			if (nodeID == masterID)
				PedigreeSearch();
			else
				PedigreeSlaveSearch();
		}
		else
			PedigreeDescribeModels();
	else
		CaseControlSearch();

#ifdef USE_MPI
	MPI::Finalize();
#endif
}


bool ModelFinder::NextTestID(int &testID) {
	bool moreWork = false;
	TESTLOCK;
	if (curTest < configuration->pTestCount)  {
		testID = curTest++;
		moreWork = true;
	}
	TESTUNLOCK;

	return moreWork;
}

#ifdef USE_MPI 

struct threadArgs {
	ModelFinder *main;
	PTestDistribution *fDist;
	PTestDistribution *oDist;
	threadArgs(ModelFinder *main, PTestDistribution *fDist, PTestDistribution *oDist) : main(main), fDist(fDist), oDist(oDist) {}
};

void *ModelFinder::ProcessWatcher(void *arg) {
	//Let's create a nice little file to write details to
	stringstream nodeLogFilename;
	
	ModelFinder *main = ((threadArgs*)arg)->main;

	nodeLogFilename<<main->configuration->reportName;
	nodeLogFilename<<".node"<<main->nodeID<<".log";
	ofstream log(nodeLogFilename.str().c_str(),ios::out|ios::trunc);


	PTestDistribution *fDist = ((threadArgs*)arg)->fDist;
	PTestDistribution *oDist = ((threadArgs*)arg)->oDist;

	int procCount = main->processCount;
	int maxModelSize = main->configuration->comboEnd;
	int activeProcesses = 0;
	int tests[2];	
	
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			if (main->NextTestID(tests[0])) {
				activeProcesses++;
				//We can't have this value being overwritten
				MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
				log<< "Test #" << tests[0] << " sent to node " <<proc<<"\n";
			}
		}
	}
	log<<"tMaster Node, "<<main->nodeID<<", initiated "<<activeProcesses<<" p-tests\n";

	//Wait for the recipients to finish and respond via round robin/
	//This thread will NEVER stop if a node can't return after being assigned
	//a legitimate test ID
	while (activeProcesses > 0) {
		float pResults[(maxModelSize + 2) * 2];

		for (int proc=0; proc<procCount; proc++) {
			pResults[proc]=0.0;
			if (proc!=main->nodeID) {		
				int idSize = 2048;
				char modelIDs[2048];
				if (MPI::COMM_WORLD.Iprobe(proc, RESULT)) {
					MPI::COMM_WORLD.Recv(tests, 2, MPI::INT, proc, RESULT);
					MPI::COMM_WORLD.Recv(pResults, (maxModelSize + 2) * 2, MPI::FLOAT, proc, RESULT);
					MPI::COMM_WORLD.Recv(modelIDs, idSize, MPI::CHAR, proc, RESULT);
					
					log<<"#"<<tests[0]<<"@Node "<<proc;
					int modelCount = main->configuration->comboEnd;
					uint resIdx = 0;

					stringstream models(modelIDs);

					for (int order = 0; order < modelCount; order++) {
						string id;
						models>>id;
						float result = pResults[resIdx++];
						if (fDist) {
							fDist->AppendTest(result, order, id.c_str(), tests[0]);
							log<<" "<<id<<"( "<<result<<"/";
						}
						result = pResults[resIdx++];
						if (oDist) {
							oDist->AppendTest(result, order, id.c_str(), tests[0]);
							log<<result<<") ";
		
						}
					}
					log<<"\n";
					if (main->NextTestID(tests[0])) {
						MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
						log<< "Test #" << tests[0] << " sent to node " <<proc<<"\n";
					}
					else
						activeProcesses--;
				}
			}
		}
	}

	//Signal each thread to go away
	tests[0] = -1;
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			log<<"Attempting to shut down node "<<proc<<"\n";
			MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
		}
	}
	log.close();
	return NULL;
	
}
#endif
/*

void ModelFinder::AppendPTest(int testID, int modelCount, float *pResult, PTestDistribution *dist, const char* modelIDs, ostream *log) {
	//Let's assume thread safety here
	//cout<<setw(6)<<testID<<setw(15)<<setprecision(3)<<pResult[0]<<setw(30)<<pResult[1]<<setw(15)<<modelCount<<"\n";
	
	for (uint modelSize = configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
		//cout<<"$  Recieved "<<pResult[modelSize+2]<<" for testID "<<testID<<" / "<<modelSize<<"\n";
		*os<<setw(10)<<pResult[pResult[
		dist->AppendTest(pResult[modelSize], modelSize, modelIDs, testID);
	}

}*/


void ModelFinder::PedigreeSlaveSearch() {
#ifdef USE_MPI
	//Let's create a nice little file to write details to
	stringstream nodeLogFilename;
	
	nodeLogFilename<<configuration->reportName<<".node"<<nodeID<<".log";
	ofstream log(nodeLogFilename.str().c_str(),ios::out|ios::trunc);

	GtFileParser *fileParser=configuration->GetInputFileParser(true);

	SnpRepository ptest;
	ptest.ParseInputFile( fileParser, NULL );

	log<<"Node "<<nodeID<<" up as slave model. Waiting for PTest data\n";
	int tests[2];

	//Blocked receive the ID we are to run. 0 Means we aren't going to do anything else
	MPI::COMM_WORLD.Recv(tests, 2, MPI::INT, masterID, TESTID);
	log<<setw(8)<<"Test #"<<setw(16)<<"Best Model ID"<<setw(8)<<"Consist."<<setw(16)<<"MDR-PDT Stat."<<setw(16)<<"Odds Ratio"<<"\n";
	while (tests[0] >= 0) {
		float pResults[(configuration->comboEnd + 2)*2];
		timer progress;

		SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);

		//This is a quick fix for random numbers. 
		srand(configuration->pTestSeed + tests[0]);
		//Acquire the evaluation suite- this time, we want randomized status
		SnpEvalSuite *pEval = configuration->GetPdtEval(true, tests[0]);

		//Populate the test repository
		ptest.LoadData( fileParser);
		//Acquire the status. Technically, this should not change, but, in case the order changed for some reason, we want to make 
		//sure we are up to date.
		CaseControlStatus stat;
		fileParser->GetStatusMask(&stat);

		pResults[0] = (float)progress.elapsed();
		progress.restart();
		//Perform search on test data
		tests[1] = ptest.Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, pEval);

		pResults[1] = (float)progress.elapsed();
		stringstream modelIDs;
	
		bestScore.Sort();
		uint idx = 0;
		//Let's put our results into the array
		for (uint modelSize= 0; modelSize < configuration->comboEnd; modelSize++) {
			if (modelSize >= configuration->comboStart-1) {
				SnpAligned* bestModel = NULL;
				ModelStatistics st = pEval->GetTopModel(modelSize, bestModel);
				st.GatherTestStatistics(&ptest, (EvalBalancedAccuracyPDT *)pEval->GetTrainer(), NULL);
				pResults[idx++] = st.GetAvgTesting();
				pResults[idx++] = st.GetOddsRatio();
				log<<setw(8)<<tests[0]<<setw(16)<<st.label<<setw(8)<<st.xvConsistency<<setw(16)<<st.GetAvgTesting()<<setw(16)<<st.GetOddsRatio()<<"\n";

				modelIDs<<st.label<<" ";
				//modelIDs<<bestScore.GetReportModel(modelSize, 0)->GetLabel()<<" ";
			}				
		}

		MPI::COMM_WORLD.Send(tests, 2, MPI::INT, masterID, RESULT);
		MPI::COMM_WORLD.Send(pResults, (configuration->comboEnd + 2)*2, MPI::FLOAT, masterID, RESULT);
		MPI::COMM_WORLD.Send(modelIDs.str().c_str(), strlen(modelIDs.str().c_str()), MPI::CHAR, masterID, RESULT);
		MPI::COMM_WORLD.Recv(tests, 1, MPI::INT, masterID, TESTID);
	}
				
	log<<"@\t\t\tNode "<<nodeID<<" shutting down\n";
	log.close();
#else
	cout<<"Slave node IDentified but this application is not MPI enabled\n";
#endif
	
}



void ModelFinder::PedigreeDescribeModels() {
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

	}
}


#include "genetics/ptestdistribution.h"

void ModelFinder::PedigreeSearch() {
	curTest = 0;


	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser(true);

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );
	//Perform some basic reporting
	configuration->GenerateReport(cout);

	if (snps.GetSnpCount() == 0) {
		cout<<"\nNo data to search over- Take a look at the pedigree report for some tips on what was found (if anything) in your pedigree data.\n";	
		cout<<"If this file is empty, your data file ("<<configuration->inputFile<<") wasn't found. If there seem to be a reasonable number of families ";
		cout<<"listed, make sure AFFECTED_VALUE matches what is found in the data.\n";
		abort();
	}
	timer progress;

#ifdef USE_MPI
	cout<<setw(45)<<right<<"MPI Master Node: "<<nodeID<<"\n";
	cout<<setw(45)<<right<<"MPI Slave Node Count: "<<processCount<<"\n";
#endif

	//Close the locus log, if it exists
	if (configuration->logLocus)
		configuration->logLocus->Close();

	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	//Create a report to store the top models
	SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);
	SnpEvalSuite *e = configuration->GetPdtEval( false, 0 );
	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)e->GetTrainer();

	//The distribution object needs to be passed to the various ptests
	PTestDistribution *fDist = e->GetFitnessDist();
	PTestDistribution *oDist = e->GetOddsRatioDist();
	bestScore.distribution = fDist;
	BasicLog *reportLog = configuration->GetReportLog();

	ostream *sigReport = reportLog->GetStream();

//Here, we don't need to do this unless we are using MPI
#ifdef USE_MPI
	pthread_t watcherThread; 
	threadArgs args(this, fDist, oDist);	
	pthread_create(&watcherThread, NULL, &ProcessWatcher, (void*)&args);
#endif

	sleep(5);

	//If some analysis is required, do it now
	if (e) {
		progress.restart();
		//Perform the search
		if (EvalBalancedAccuracyPDT::Verbose)
			cout<<setw(16)<<right<<"Model"<<setw(8)<<"   T Statistic"<<endl;
		else
			cout<<"  Searching\n";

		uint modelCount = snps.Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, e);

		evalPdt->DisplayAnalyses(true);
		cout<<"\n\nPedigree Search of "<<e->GetSnpCount()<<" : "<<modelCount<<" completed in "<<progress.elapsed()<<" seconds\n";


		cout<<"\n\nPerforming "<<configuration->pTestCount<<" Permutation Tests:\n";
		
		//GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();
		cout<<setw(6)<<"Test #"<<setw(15)<<"Load Time (s)"<<setw(30)<<"Execution Time(s)"<<setw(15)<<"Models Seen"<<endl;

		int i =0; 
		while (NextTestID(i)) {
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
			
			for (uint modelSize= configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
				SnpAligned* bestModel = NULL;
				ModelStatistics st = pEval->GetTopModel(modelSize, bestModel);
				st.GatherTestStatistics(&ptest, (EvalBalancedAccuracyPDT *)pEval->GetTrainer(), NULL);
				fDist->AppendTest(st.GetAvgTesting(), modelSize, st.GetLabel().c_str(), i);
				oDist->AppendTest(st.GetOddsRatio(), modelSize, st.GetLabel().c_str(), i);
			}
			//pEval->UpdateDistribution();
			cout<<setw(30)<<progress.elapsed()<<setw(15)<<modelCount<<"\n";
			delete pEval;
		}

#ifdef USE_MPI
	cout<<"$  Waiting for the watcher thread to join\n";
	pthread_join(watcherThread, NULL); 
	cout<<"$ Completed. We can now report the final results\n";
	//((Genetics::Distribution::PTestDistribution*)fDist)->Report(&cout);
	//fDist->Report(&cout);
	//oDist->Report(&cout);
#endif


		//Set up the distribution and sort the results
	//	e->BuildDistribution();

		//Report on the distribution
		if (fDist)
			fDist->DumpDistribution(configuration->GetDistributionLog()->GetStream());
		if (oDist) 
			oDist->DumpDistribution(configuration->GetDistributionLog()->GetStream());
	
		float minPVal = configuration->pValueThreshold;	
		
		bestScore.Sort();
		if (sigReport == NULL) {
			cout<<"\n\n";
			e->ReportResults(&cout);
		}
		else {
			if (fDist && configuration->crossValidationCount == 1 ){
				//This is the significant model report
				for (uint modelSize = configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
					uint padding = 4 * modelSize + 4;
					uint  reportSize=bestScore.GetTotalEntryCount(modelSize, 0); //configuration->crossValidationCount);
					*sigReport<<"\n\nDistribution Report "<<modelSize + 1<<" SNP(s) per model"<<endl;
					*sigReport<<"\t* Only models with a p-value < "<<setprecision(3)<<minPVal<<" are reported. "<<endl;
					*sigReport<<setw(padding)<<""<<setw(12)<<"T    "<<setw(16)<<"Statistical"<<endl;
					*sigReport<<setw(padding)<<"Model"<<setw(12)<<"Statistic"<<setw(16)<<"Significance"<<endl;
					
					float pValue=1.0;

					uint foldID = 0;
				
					SnpReposTxtSorted::ReportModel *model = bestScore.GetReportModel(modelSize,foldID, 0);
					if (model) 
						pValue=fDist->GetPValue(model->score, modelSize);
		
					bool significantFinds=false;
					for (uint idx=1; model && idx < reportSize; idx++) {
						if (pValue <= minPVal) {
							significantFinds=true;
							*sigReport<<setw(padding)<<model->GetLabel();
							*sigReport<<setw(12)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3)<<model->score;
							*sigReport<<"       p < "<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(3)<<pValue<<endl;
						}
						model = bestScore.GetReportModel(modelSize, foldID, idx);
						if (model)
							pValue=fDist->GetPValue(model->score, modelSize);
						else
							pValue = 0.0;
						
					}
					if (!significantFinds)
						*sigReport<<setw(12)<<"None found"<<endl;
				}
			}				
			//Report on the overal resutlts
			*sigReport<<"\n\n";

			evalPdt->ReportResults(sigReport, configuration->reportModelCount, &snps);
		}
		delete e;

	}
}
void ModelFinder::CaseControlSearch() {
	//Acquire the 
	GtFileParser *fileParser=configuration->GetInputFileParser(true);
	if (configuration->logLocus)
		configuration->logLocus->WriteHeader();

	timer progress;
	snps.ParseInputFile( fileParser, configuration->logLocus );
	cout<<"\n\nInput Parse completed in "<<progress.elapsed()<<" seconds\n";

	configuration->GenerateReport(cout);
	if (configuration->logLocus)
		configuration->logLocus->Close();

	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);
	cout<<"Affected:  "<<stat.affected<<"\n";
	cout<<"Unaffected "<<stat.unaffected<<"\n";
	cout<<"Total      "<<stat.total<<"\n";
	
	cout<<"\n-----------------------------------\n";
	cout<<  "- Overall Status:                 -\n";
	cout<<  "- Affected   :"<<stat.affected.count()<<"\n";
	cout<<  "- Unaffected :"<<stat.unaffected.count()<<"\n";
	cout<<  "- Total      :"<<stat.total.count()<<"\n";
	
	SnpRepository *notPerfect = NULL;

	//Create a report to store the top models
	SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);
	SnpEvalSuite *e = configuration->GetEval( stat, fileParser);


	//If some analysis is required, do it now
	if (e) {
		if (configuration->DoPerfectSearch()) {
			//CaseControlStatus completeStatus = stat.CombinedStatus();
			IsPerfect perfect(stat);
		
			//Create the report to store perfect models
			SnpReposTxtSorted isPerfect(0, 1);
			configuration->CloseInputParser();
		
			//This is where we will filter the imperfect models
			notPerfect=new SnpRepository();
			notPerfect->InitRepository(configuration->affectedCount +configuration->unaffectedCount, 0);
			progress.restart();	

			snps.Evaluate( &isPerfect, notPerfect, &perfect);
			cout<<"Perfect Search of "<<perfect.GetSnpCount()<<" models completed in "<<progress.elapsed()<<" seconds\n";
			//configuration->perfReport->ReportOnRepository( "Perfect 1-SNP Models", 1, &isPerfect);	
		}
		else
			notPerfect = &snps;
		progress.restart();



		notPerfect->Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, e);

		cout<<"\n\nMD Search of "<<e->GetSnpCount()<<" completed in "<<progress.elapsed()<<" seconds\n";

		e->BuildDistribution();

		if (e->GetFitnessDist())
			e->GetFitnessDist()->DumpDistribution(configuration->GetDistributionLog()->GetStream());

		e->ReportResults(&cout);

		//only delete it if it's a new thing
		if (configuration->DoPerfectSearch())
			delete notPerfect;
	} //Else, we have nothing really to do to change the repository
	else {
		AnySnp any;
		configuration->CloseInputParser();
		SnpRepository *results=new SnpRepository();
		results->InitRepository(configuration->affectedCount +configuration->unaffectedCount, 0);
		snps.Evaluate( results, &any);		
	}
}




}
