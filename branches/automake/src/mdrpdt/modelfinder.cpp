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
//#include "snprepostxtsorted.h"
#include <boost/timer.hpp>





#ifndef TEST_APP
int main(int argc, char **argv) {
	MDR::ModelFinder app;

	if (app.ParseCmdLine(argc, argv) && app.PerformSanityCheck())
		app.Start();
	return 0;
}
	
#endif

namespace MDR {

using namespace boost;

ModelFinder::ModelFinder() : doAnalyze(true) {
#ifdef USE_MPI
	appname="parallel mdr-pdt";
#else
	appname="mdr-pdt";
#endif
	appfunction="Multifactor Dimensionality Reduction - Pedigree Disequilibrium Test";
	authors="Marylyn Ritchie & Eric Torstenson\nPlease forward any comments or errors to: mdr-pdt@chgr.mc.vanderbilt.edu";
	major=2;
	minor=5;
	bugFixes=5;
	configuration = EseConfiguration::Instance();

	rand = &Utility::Random::globalGenerator;
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
	ss<<" * Set the number of cross validation folds (we recommend 1, 5 or 10. 1 Is normal execution)\n";
	ss<<"CROSSVALINTERVAL       1\n";
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
		ss<<"\n\n******************* Mendelian Error Detection\n";
		ss<<"* Depending on the data provided, mendelian errors can be detected and filtered out according to user preference\n";
		ss<<"* 3 Variables help control the action taken when errors are discovered. It should be noted that parental data is \n";		
		ss<<"* present before errors can be detected.\n";
		ss<<"* \n";
		ss<<"* MENDELIAN_ERROR_LEVEL level\n";
		ss<<"* This setting gives the user 3 levels of control.\n";
		ss<<"*   1 - This tells the application to report all errors but do nothing to the data before analyzing. It should be \n";
		ss<<"*       noted that your power will decrease with the presence of mendelian errors\n";
		ss<<"*   2 - This instructs the application to zero out all individuals in the pedigree at that locus. This is using the\n";
		ss<<"*       idea that the user doesn't really know where the error occurred. \n";
		ss<<"*   3 - In addition to stripping the data from all members of the pedigree at the specified locus, the application does the following:\n";
		ss<<"*       *   Evaluate error counts per pedigree and remove the whole pedigree if the error count exceeds the pedigree threshold. \n";
		ss<<"*       *   Evaluate errors at each locus and zero out all individuals at loci where the number of errors exceed the locus threshold.\n";
		ss<<"MENDELIAN_ERROR_LEVEL 1\n";
		ss<<"*   Indicates the threshold for zeroing out all data at a given locus (percentage of all pedigrees where 1 or more errors occurred\n";
		ss<<"*   at the locus\n";
		ss<<"MENDELIAN_LOCUS_THRESHOLD 25\n";
		ss<<"*   Indicate the percentage of total genotypes with errors are required before the pedigree is dropped from consideration\n";
		ss<<"MENDELIAN_PEDIGREE_THRESHOLD 25\n";
		ss<<"*   Indicate whether or not a version of the 'cleaned' data is written to file\n";
		ss<<"WRITE_CLEAN_DATAFILE Off\n";


		ss<<"\n\n/****************** Input format\n";
		ss<<" * PEDIGREE   - For PDT anayses there is only one format supported\n";
		ss<<"INPUTFORMAT            PEDIGREE\n";
		ss<<" *  The name of the input file where your SNP data is to be found. This file must be space delimited\n";
		string filename = "YourFilename.ped";
		if (argc>2)
			filename=argv[2];
		ss<<"INPUTFILE              "<<filename<<"\n";
		specialExt<<"EXT_PEDIGREE          pedigree      //Pedigree related errors and notes are logged here\n";
		specialExt<<"EXT_CLEANED_DATA      No            //If data is cleaned, do we output it as a new file?\n";
		ss<<" *  There is only one Pedigree analysis currently available\n";
		
		ss<<" *  You can exclude certain pedigrees from analysis. Buy removing the asterisk from the line below\n";
		ss<<" *  will tell the application to exclude pedigrees 13, 21 and 1003\n";
		ss<<"*PEDIGREE_EXCLUSIONS 13 21 1003\n";

		ss<<"ANALYSISSTYLE          PDT\n";
		ss<<" *  Turn On/Off Triad expansion\n";
		ss<<"EXPAND_TRIOS           On\n";
		ss<<" *  Exapand all affected (when parental data exists). This doesn't just apply to trios\n";
		ss<<"EXPAND_ALL             Off\n";
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
		ss<<" *  Writes the arrangement of the pedigree assignments (for cross validation folding) to the pedigree report\n";
		ss<<"VERBOSE_FOLDING        No\n";
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
#ifdef USE_MPI
	cout<<"usage: mdr-pdtp <configuration file> [optional model IDs]\n";
#else
	cout<<"usage: mdr-pdts <configuration file> [optional model IDs]\n";
#endif
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

	MpiPrep();

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

	MpiClose();
}




#ifdef USE_MPI 

struct mdrPdtThreadArgs {
	ModelFinder *main;
	PTestDistribution *fDist;
	PTestDistribution *oDist;
	PTestDistribution *pDist;
	PTestDistribution *lrDist;
	mdrPdtThreadArgs(ModelFinder *main, PTestDistribution *fDist, PTestDistribution *oDist, PTestDistribution *pDist, PTestDistribution *lrDist) 
			: main(main), fDist(fDist), oDist(oDist), pDist(pDist), lrDist(lrDist) {
		cout<<"mdrPdtThreadArgs Constructed\n";
	}

	~mdrPdtThreadArgs() {
		cout<<"mdrPdtThreadArgs Destructed\n";
	}
};

void *ModelFinder::ProcessWatcher(void *arg) {
	//Let's create a nice little file to write details to
	stringstream nodeLogFilename;
	
	ModelFinder *main = ((mdrPdtThreadArgs*)arg)->main;

	nodeLogFilename<<main->configuration->reportName;
	nodeLogFilename<<".node"<<main->nodeID<<".log";
	ofstream log(nodeLogFilename.str().c_str(),ios::out|ios::trunc);


	PTestDistribution *fDist = ((mdrPdtThreadArgs*)arg)->fDist;
	PTestDistribution *oDist = ((mdrPdtThreadArgs*)arg)->oDist;
	PTestDistribution *pDist = ((mdrPdtThreadArgs*)arg)->pDist;
	PTestDistribution *lrDist = ((mdrPdtThreadArgs*)arg)->lrDist;

	int procCount = main->processCount;
	int maxModelSize = main->configuration->comboEnd;
	int activeProcesses = 0;
	size_t tests[2];	
	int masterID = main->masterID;
	int nodeID = main->nodeID;
	int buffSize = (maxModelSize + 3) * 2;
	
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != masterID) {
			if (main->NextTestID(tests[0])) {
				activeProcesses++;
				MPILOCK;
				//We can't have this value being overwritten
				MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
				MPIUNLOCK;
				log<< "Test #" << tests[0] << " sent to node " <<proc<<"\n";
				log.flush();
			}
		}
	}
	log<<"Master Node, "<<nodeID<<", initiated "<<activeProcesses<<" p-tests\n";
	log.flush();
	//Wait for the recipients to finish and respond via round robin/
	//This thread will NEVER stop if a node can't return after being assigned
	//a legitimate test ID
	while (activeProcesses > 0) {
		sleep(3);

		float pResults[buffSize];

		for (int proc=0; proc<procCount; proc++) {
			pResults[0]=0.0;
			log.flush();
			if (proc!=nodeID) {		
				log.flush();
				int idSize = 2048;
				char modelIDs[2048];  
				log.flush();
				MPILOCK;
				if (MPI::COMM_WORLD.Iprobe(proc, RESULT)) {
					//assert(proc != masterID);
					MPI::COMM_WORLD.Recv(tests, 2, MPI::INT, proc, RESULT);
					MPI::COMM_WORLD.Recv(pResults, buffSize, MPI::FLOAT, proc, RESULT);
					MPI::COMM_WORLD.Recv(modelIDs, idSize, MPI::CHAR, proc, RESULT);
					
					log<<"#"<<tests[0]<<"@Node "<<proc;
					log.flush();

					int resIdx = 0;

					stringstream models(modelIDs);

					for (int order = 0; order < maxModelSize; order++) {
						string id;
						models>>id;
						if (resIdx + 1 > buffSize) {
							log<<"  -- Trying to work outside buffer limits "<<resIdx<<", "<<resIdx+1<<" : "<<buffSize<<"\n";
							log.flush();
						}
						float result = pResults[resIdx++];
						if (fDist) {
							log<<" "<<id<<"( "<<result;
							log.flush();
							fDist->AppendTest(result, order, id.c_str(), tests[0] - 1);
						}
						result = pResults[resIdx++];
						if (oDist) {
							log<<"/"<<result;
							log.flush();		
							oDist->AppendTest(result, order, id.c_str(), tests[0] - 1);
						}
						result = pResults[resIdx++];
						if (pDist) {
							log<<"/"<<result;
							log.flush();
							pDist->AppendTest(result, order, id.c_str(), tests[0] - 1);
						}
						
						log<<") ";
						log.flush();
					}
					log<<"\n";
					log.flush();
					if (main->NextTestID(tests[0])) {
						MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
						log<< "Test #" << tests[0] << " sent to node " <<proc<<"\n";
						log.flush();
					}
					else
						activeProcesses--;
				}
				MPIUNLOCK;
			}
		}
	}

	MPILOCK;
	//Signal each thread to go away
	tests[0] = -1;
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			log<<"Attempting to shut down node "<<proc<<"\n";
			MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
		}
	}
	MPIUNLOCK;
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

	fileParser=configuration->GetInputFileParser(true);

	//We need to initiate the load- but, we don't want to keep the set after each run
	SnpRepository loader;					
	loader.ParseInputFile( fileParser, NULL );


	log<<"Node "<<nodeID<<" up as slave model. Waiting for PTest data\n";
	int tests[2];

	MPILOCK;
	//Blocked receive the ID we are to run. 0 Means we aren't going to do anything else
	MPI::COMM_WORLD.Recv(tests, 2, MPI::INT, masterID, TESTID);
	MPIUNLOCK;

	log<<setw(8)<<"Test #"<<setw(16)<<"Best Model ID"<<setw(8)<<"Cons.";
	log<<setw(16)<<"MDR-PDT Stat."<<setw(16)<<"Odds Ratio"<<"\n";
	log.flush();

	int buffSize = (configuration->comboEnd + 3)*2;
	stringstream modelIDs;

	while (tests[0] >= 0) {
		float pResults[buffSize];
		vector<ModelStatistics> bestModels;
		pResults[0] = 0.0;				//I can't really worry about the load time- too many trivial values already
		pResults[1] = RunPTest(tests[0], tests[1], bestModels);
		
		uint idx = 0;
		int stIdx = 0;
		
		//Let's put our results into the array. We will be iterating over model orders
		for (uint modelSize= 0; modelSize < configuration->comboEnd; modelSize++) {
			if (modelSize >= configuration->comboStart-1) {
				ModelStatistics &st = bestModels[stIdx++];
				pResults[idx++] = st.GetAvgTesting();
				pResults[idx++] = st.GetOddsRatio();
				pResults[idx++] = st.GetAvgPredictionError();
				log<<setw(8)<<tests[0]<<setw(16)<<st.label<<setw(8)<<st.xvConsistency;
				log<<setw(16)<<st.GetAvgTesting()<<setw(16)<<st.GetOddsRatio()<<"\n";
				log.flush();
				modelIDs<<st.label<<" ";
				//modelIDs<<bestScore.GetReportModel(modelSize, 0)->GetLabel()<<" ";
			}				
		}		
		MPILOCK;
		MPI::COMM_WORLD.Send(tests, 2, MPI::INT, masterID, RESULT);
		MPI::COMM_WORLD.Send(pResults, (configuration->comboEnd + 2)*2, MPI::FLOAT, masterID, RESULT);
		MPI::COMM_WORLD.Send(modelIDs.str().c_str(), strlen(modelIDs.str().c_str())+1, MPI::CHAR, masterID, RESULT);
		MPI::COMM_WORLD.Recv(tests, 2, MPI::INT, masterID, TESTID);
		MPIUNLOCK;

	}

				
	log<<"@\t\t\tNode "<<nodeID<<" shutting down\n";
	log.close();
#else
	cout<<"Slave node IDentified but this application is not MPI enabled\n";
#endif
	
}



void ModelFinder::PedigreeDescribeModels() {
	//Acquire the input parser from the configuration
	Utility::Random::globalGenerator.Seed(configuration->pTestSeed);

	fileParser=configuration->GetInputFileParser(false);
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();

	snps.ParseInputFile( fileParser, NULL );
	configuration->ValidateSnpCount(snps.GetSnpCount());


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
	curTest = 1;
	Utility::Random::globalGenerator.Seed(configuration->pTestSeed);

	//Acquire the input parser from the configuration
	fileParser=configuration->GetInputFileParser(true);

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );

	configuration->ValidateSnpCount(snps.GetSnpCount());

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
	PTestDistribution *pDist = e->GetPEDist();
	bestScore.distribution = fDist;

	BasicLog *reportLog = configuration->GetReportLog();

	ostream *sigReport = reportLog->GetStream();

//Here, we don't need to do this unless we are using MPI
#ifdef USE_MPI
	pthread_t watcherThread; 
	mdrPdtThreadArgs args(this, fDist, oDist, pDist, NULL);	
	pthread_create(&watcherThread, NULL, &ProcessWatcher, (void*)&args);
#endif

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
		
		cout<<setw(6)<<"Test #"<<setw(15)<<"Load Time (s)"<<setw(30)<<"Execution Time(s)"<<setw(15)<<"Models Seen"<<endl;


		//Let's bump the test up so we can know we are working on PTests
		size_t i =0; 
		while (NextTestID(i)) {
			int modelCount;
			vector<ModelStatistics> topModels;
			float searchTime = RunPTest(i, modelCount, topModels);

			uint stIdx = 0;
			cout<<"PTest "<<i<<" ";
			for (uint modelSize= 0; modelSize < configuration->comboEnd; modelSize++) {
				if (modelSize >= configuration->comboStart-1) {
					//If this fails the RunPTest probably failed
					assert(topModels.size() > stIdx);

					ModelStatistics &st = topModels[stIdx++];
					st.GatherTestStatistics(&snps, evalPdt, NULL);

					cout<<st.GetLabel()<<" ";
					fDist->AppendTest(st.GetAvgTesting(), modelSize, st.GetLabel().c_str(), i - 1);
					oDist->AppendTest(st.GetOddsRatio(), modelSize, st.GetLabel().c_str(), i - 1);
					pDist->AppendTest(st.GetAvgPredictionError(), modelSize, st.GetLabel().c_str(), i - 1);
				}
			}			
//			cout<<"\n";
			cout<<setw(30)<<searchTime<<setw(15)<<modelCount<<"\n";

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
		if (pDist)
			pDist->DumpDistribution(configuration->GetDistributionLog()->GetStream());	

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
					uint  reportSize=bestScore.GetTotalEntryCount(0, modelSize); //configuration->crossValidationCount);
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
