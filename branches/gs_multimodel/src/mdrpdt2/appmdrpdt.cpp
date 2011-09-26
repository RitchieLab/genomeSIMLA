//
// C++ Implementation: appmdrpdt
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "appmdrpdt.h"
#include "tstatistic.h"
#include "utility/exception.h"
#include "pedigreecleaner.h"
#include "pedigreeexclusion.h"
#include "pedigreestatistics.h"
#include <pthread.h>
#include "missingdataevaluation.h"
#include "timestamp.h"

using namespace MdrPDT::Validator;

namespace MdrPDT {
#ifdef USE_MPI
pthread_mutex_t AppMdrPDT::input_lock 	= PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t AppMdrPDT::output_lock 	= PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t AppMdrPDT::mpi_lock		= PTHREAD_MUTEX_INITIALIZER;
#endif
pthread_mutex_t ptestLock = PTHREAD_MUTEX_INITIALIZER;

AppMdrPDT::AppMdrPDT() 
		: seed(1397), xvCount(5), eval(NULL), winningScore(0.0), 
		ptestShortCircuit(0), ptestsRequested(0), ptestExceptions(0), 
		threadCount(2), nodeID(0), masterID(0), processCount(0) {	
}


AppMdrPDT::~AppMdrPDT() {	
	MpiClose();
	if (eval)
		delete eval;
}

void AppMdrPDT::PrintBanner() {
	cout<<"mdrpdt "<<APPMAJOR<<"."<<APPMINOR<<"."<<APPBUGFIX<<" ("<<BUILD_NUMBER<<") "<<BUILD_TYPE<<"  "<<BUILD_DATE<<"\n";
	cout<<"Multifactor Dimension Reduction - Pedigree Disequilibrium Test\n";
#ifdef USE_MPI
	cout<<"* This application is compiled to run on parallel computing systems using MPI\n";
#else
	cout<<"* (serial)\n";
#endif
	cout<<"Marylyn Ritchie, Eden Martin, Todd Edwards and Eric Torstenson\nPlease forward any comments or errors to mdr-pdt@chgr.mc.vanderbilt.edu\n";
	
}

void AppMdrPDT::PrintHelp() {
	PrintBanner();
#ifdef USE_MPI
	cout<<"usage: mdr-pdtp <configuration file> [optional model IDs]\n";
#else
	cout<<"usage: mdrpdt <configuration file> [optional model IDs]\n";
#endif
	cout<<"\nThis search tool can be used to perform analyses on SNP data or to explore a single model\n";
	cout<<"* To execute analysis, simply run the command with the configuration file only\n";
	cout<<"* To explore a specific model, specify the configuration file and one or more models separated by space\n";
	cout<<"* To dump an example configuration file, use the following format\n";
	cout<<"*\t\tmdrpdt PDT [ped file] [min model size] [max model size] [xv count] [ptest count] [random seed]\n";
	cout<<"*\t\tEach parameter after PDT is optional, however, they must be in the order specified above.\n"; 
}


bool AppMdrPDT::ParseCmdLine(int argc, char **argv) {
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif

	if (argc < 2) {
		PrintHelp();
		return false;
	}

	//PDT [cfg] [start] [stop] [xv] [ptest] [seed]
	if (strcmp(argv[1], "PDT") == 0) {
		cfg.Load();
		if (argc > 2) {
			//First argument is the pedigree file
			cfg.appConfiguration->SetValue("INPUTFILE", string(argv[2]));
		}
		if (argc > 4) {
			cfg.appConfiguration->SetValue("COMBO_START", atoi(argv[3]));
			cfg.appConfiguration->SetValue("COMBO_END", atoi(argv[4]));
		}
		if (argc > 5) {
			cfg.appConfiguration->SetValue("CROSSVALINTERVAL", atoi(argv[5]));
		}
		if (argc > 6) {
			cfg.appConfiguration->SetValue("PTEST_COUNT", atoi(argv[6]));
		}
		if (argc > 7){ 
			cfg.appConfiguration->SetValue("PTEST_SEED", atoi(argv[7]));
		}
		cfg.appConfiguration->WriteConfiguration(cout);
		exit(0);
	}
	LoadConfiguration(argv[1]);
	cfg.ReportConfiguration(cout);

	if (argc > 2) {
		targetModel = argv[2];
	}
	return true;
}
ApplicationConfiguration *AppMdrPDT::LoadConfiguration(const char *cfgFilename) {
	ApplicationConfiguration *appConfig = cfg.Load(cfgFilename);
	appConfig->ExecuteConfiguration();
	seed = appConfig->GetInteger("PTEST_SEED");
	xvCount = appConfig->GetInteger("CROSSVALINTERVAL");
	threadCount = appConfig->GetInteger("THREAD_COUNT");

	MpiPrep();
	return appConfig;
}

string AppMdrPDT::GetFilename(const char *extension) {
	return cfg.GetFilename() + "." + string(extension);
}

void AppMdrPDT::LoadData(const char *dataset) {
	pedigreeReport.open(GetFilename(cfg.appConfiguration->GetLine("EXT_PEDIGREE").c_str()).c_str());
	
	//Family information. This is used to generate 
	try {
		repo.Load(dataset, pedigreeReport);
		cout<<setw(45)<<"Configuration File Name: "<<cfg.GetFilename()<<"\n";
		cout<<setw(45)<<"Dataset File Name: "<<dataset<<"\n";
		cout<<setw(45)<<"Total Number of Loci in Dataset: "<<repo.GetLocusCount()<<"\n";
	} catch (Utility::Exception::FileNotFound &e) {
		cerr<<"Unable to load dataset, "<<dataset<<".\n";
	}

	if (Individual::UseMerlin) {
		string datFilename = cfg.appConfiguration->GetLine("MERLIN_DAT");
		if (datFilename.length() > 0) {
			repo.LoadDat(cfg.appConfiguration->GetLine("MERLIN_DAT").c_str());
			cout<<setw(45)<<"DAT File: "<<cfg.appConfiguration->GetLine("MERLIN_DAT")<<"\n";
		}
	}

	//First, let's exclude pedigrees
	vector<string> pedigrees;
	if (cfg.appConfiguration->GetLines("EXCLUDE_PEDIGREES", pedigrees)) {
		PedigreeExclusion pedEx(pedigrees);
		pedEx.EvaluateRepository(&repo);
		if (pedEx.ExclusionCount() > 0) {
			cout<<setw(45)<<"Excluded Pedigrees: "<<pedEx.ExclusionCount()<<"\n";
			cout<<setw(45)<<"Exclusions: "<<pedEx.ExcludedPedigrees()<<"\n";
		}
		else
			cout<<setw(45)<<"Excluded Pedigrees: "<<"None"<<"\n";
	}
	else
		cout<<setw(45)<<"Excluded Pedigrees: "<<"None"<<"\n";

	//Next, we run the cleaner on the data
 	Validator::PedigreeCleaner cleaner(cfg.appConfiguration->GetInteger("MENDELIAN_ERROR_LEVEL"), cfg.appConfiguration->GetInteger("MENDELIAN_PEDIGREE_THRESHOLD"));
	cleaner.Evaluate(&repo, pedigreeReport);

	int totalDSPs = repo.PostLoad(pedigreeReport);
	//Finally, we will get some data about what is left and report it
	Validator::PedigreeStatistics stats(cout);
	stats.EvaluateRepository(&repo);
	cout<<setw(45)<<"Participating DSPs: "<<totalDSPs<<"\n";
	
	

	if (cfg.appConfiguration->GetBoolean("WRITE_CLEAN_DATAFILE")) 
		repo.Write(GetFilename(cfg.appConfiguration->GetLine("EXT_CLEANED_DATA").c_str()).c_str());

}

Evaluation::EvaluationMethod *AppMdrPDT::GetEvaluationMethod() {
	return eval;
}
void AppMdrPDT::SetSeed(int seed) { 
	this->seed = seed; 
	cfg.appConfiguration->SetValue("PTEST_SEED", seed);
}
void AppMdrPDT::EvaluatePTest(int ptest, Evaluation::TFinalReport& report) {
	Utility::Random rnd(seed+ptest);
	//Initialize the data for real analysis, not permuted
	repo.InitializeData(data, rnd, xvCount, true);
	
	eval = new Evaluation::TStatistic(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());

	eval->BasicEval(&data, report);
}


void AppMdrPDT::SetBestModel(float bestMOR, const char *bestModel, int ptestCount) {
	this->winningScore = bestMOR;
	this->winningID = bestModel;
	float threshold = cfg.appConfiguration->GetDouble("PTEST_SHORTCIRCUIT");
	if (threshold > 0.0) {
		this->ptestShortCircuit = (int)(threshold * ptestCount);
		cout<<"\n\nPTests will short circuit if "<<ptestShortCircuit<<" tests exceed value: "<<bestMOR<<"\n";
	}
}

void AppMdrPDT::EvaluateDataset(Evaluation::TFinalReport& report)   {
	//The report which keeps all of our results from the real run
	report = Evaluation::TFinalReport(
			Evaluation::EvaluationMethod::maxModelSize, 
			xvCount, 
			cfg.appConfiguration->GetInteger("REPORTMODELCOUNT"));

	Utility::Random rnd(cfg.appConfiguration->GetInteger("PTEST_SEED"));

	//Initialize the data for real analysis, not permuted
	repo.InitializeData(data, rnd, xvCount, false);
	
	eval = new Evaluation::TStatistic(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());
	Evaluation::MissingDataEvaluation mde(cfg.appConfiguration->GetDouble("MISSING_THRESHOLD"));

	//Work through specific model. There are a few things that we don't worry about
	if (targetModel.length() > 0) {
		eval->EvaluateModel(&data, targetModel.c_str());
		int lociInConsideration = mde.Analyze(cout, data);
		cout<<setw(45)<<"Loci For Analysis: "<<lociInConsideration<<"\n";
	}
	//Work through a complete analysis
	else {
		vector<string> exclusions;
		data.InitExclusionList(exclusions);
		int lociInConsideration = mde.Analyze(cout, data);
		cout<<setw(45)<<"Loci For Analysis: "<<lociInConsideration<<"\n";
		if (lociInConsideration < Evaluation::EvaluationMethod::maxModelSize) {
			cerr<<"After removing bad loci and those from the exclusion list from the list of SNPs, there are not enough to perform the desired analysis. Please check that your data is properly configured\n";
			exit(1);
		}
		eval->BasicEval(&data, report);
	}

}





void AppMdrPDT::ResetPTestCounts(int ptestCount) {
	this->ptestCount = ptestCount;
	ptestExceptions=0;
	ptestsRequested=0;
}

struct thread_params {
	Distribution::PTestDistribution *dist;
	AppMdrPDT *app;

	thread_params(Distribution::PTestDistribution* dist, AppMdrPDT* app) : dist(dist), app(app) { }
};

void *AppMdrPDT::RunPTests(void *args) {
	thread_params *p=(thread_params*)args;
	p->app->RunPTests(p->dist);
	delete p;
	return NULL;
}




#ifndef USE_MPI /***************** non-MPI Versions ***********************/
int AppMdrPDT::GetNextPTest() {
	pthread_mutex_lock(&ptestLock);
	int returnVal = -1;
	if (ptestsRequested < ptestCount)
		returnVal = ++ptestsRequested;
	pthread_mutex_unlock(&ptestLock);
	return returnVal;
}
Distribution::PTestDistribution *AppMdrPDT::RunPTests(int ptestCount) {
	Distribution::PTestDistribution *morDist = new Distribution::OmnibusDistribution(ptestCount);
	int maxThreadCount=threadCount-1;
	pthread_t threads[maxThreadCount]; 

	ResetPTestCounts(ptestCount);
	
	//Kick off the threads (leaving 1 for the main process)
	for (int i=0; i<maxThreadCount; i++) 
		pthread_create(&threads[i], NULL, RunPTests, (void*)new thread_params(morDist, this));
	RunPTests(morDist);
	for (int i=0; i<maxThreadCount; i++) 
		pthread_join(threads[i], NULL);
	return morDist;
}

void AppMdrPDT::RunPTests(Distribution::PTestDistribution* disto) {
	pthread_mutex_lock(&ptestLock);
	Evaluation::TStatistic eval(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());
	pthread_mutex_unlock(&ptestLock);

	int testNumber = GetNextPTest();
	while (testNumber > 0) {
		Utility::Random rnd(seed+testNumber);
		GenotypeRepository test;

		//Create the test data partition
		pthread_mutex_lock(&ptestLock);
		repo.InitializeData(test, rnd, xvCount, true);
		pthread_mutex_unlock(&ptestLock);

		vector<string> exclusions;
		if (cfg.appConfiguration->GetLines("EXCLUDE_LOCUS", exclusions))
			test.InitExclusionList(exclusions);
		Evaluation::TStatistic testEval(xvCount, repo.GetPedigreeCount(), test.GetIndividualCount());
		Evaluation::TFinalReport ptestWinner(Evaluation::EvaluationMethod::maxModelSize, xvCount, 1);

		testEval.BasicEval(&test, ptestWinner);

		//We use an omnibus distribution, so we need to get the best of all N orders
		string bestID = "";
		float bestScore = 0.0;
		for (int o=Evaluation::EvaluationMethod::minModelSize-1; o<Evaluation::EvaluationMethod::maxModelSize;o++) {
			string modelID;
			float curScore = ptestWinner.GetBestMOR(o, modelID);
			if (curScore>bestScore) {
				bestScore=curScore;
				bestID = modelID;
			}
		}
			
		pthread_mutex_lock(&ptestLock);
		disto->AppendTest(1, testNumber-1, bestID.c_str(), bestScore);
		pthread_mutex_unlock(&ptestLock);


		std::cout<<".";std::cout.flush();
		if (ptestShortCircuit > 0) {
			if (bestScore>winningScore) 
				ptestExceptions++;
			if (ptestExceptions>ptestShortCircuit) {
				cerr<<"\n\nToo many permutations have outperformed the top model (Model: "<<
					winningID<<" MOR="<<winningScore<<".) Halting execution prematurely\n";
				ptestsRequested=ptestCount;
			}
		}
		testNumber = GetNextPTest();
	}
}

void AppMdrPDT::MpiPrep() {

}

void AppMdrPDT::MpiClose() {

}
#else	/************************** MPI Versions **************************/
int AppMdrPDT::GetNextPTest() {
	/**
	 * We have two ways this works.
	 * 1)	Master Node will run normally
	 * 2)	Slave Nodes will actually query Master Node for number and return
	 */
	int returnVal = -1;
	if (nodeID == masterID) {
		pthread_mutex_lock(&ptestLock);
		if (ptestsRequested < ptestCount)
			returnVal = ++ptestsRequested;
		pthread_mutex_unlock(&ptestLock);
	}
	else {
		MPILOCK;
		MPI::COMM_WORLD.Recv(&returnVal, 1, MPI::INT, masterID, MPI_TESTID);
		mpiReport<<">"<<returnVal<<"\n"; mpiReport.flush();
		MPIUNLOCK;
	}
	return returnVal;
}
struct threadDetails {
	Distribution::PTestDistribution* dist;
	AppMdrPDT *main;
	
	threadDetails(Distribution::PTestDistribution* dist, AppMdrPDT *main) : 
		dist(dist), main(main) { }
};
Distribution::PTestDistribution *AppMdrPDT::RunPTests(int ptestCount) {
	Distribution::PTestDistribution *morDist = new Distribution::OmnibusDistribution(ptestCount);
	pthread_t watcherThread;
	if (masterID == nodeID) {
		ResetPTestCounts(ptestCount);
			
		//Start the watcher in it's own thread. This will take care of all of the slave nodes
		threadDetails details(morDist, this);
		pthread_create(&watcherThread, NULL, &ProcessWatcher, (void*)&details);
	}
	RunPTests(morDist);
	if (masterID == nodeID)
		pthread_join(watcherThread, NULL);
	return morDist;
}


	

//This sends out each test to remote nodes and receives/tabulates the results
void *AppMdrPDT::ProcessWatcher(void *arg) {
	AppMdrPDT *app = ((threadDetails*)arg)->main;
	Distribution::PTestDistribution* dist = ((threadDetails*)arg)->dist;
	int procCount = app->processCount;
	//EvaluationMethod::maxModelSize;
	int activeProcesses = 0;

	int pTest = 0;				/** Used to store the test offset to be done/reported on */
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Avoid sending ourselves a message.
		if (proc != app->nodeID) {
			pTest = app->GetNextPTest();
			if (pTest >= 0) {
				activeProcesses++;
				MPILOCK;
				MPI::COMM_WORLD.Send(&pTest, 1, MPI::INT, proc, MPI_TESTID);
				app->mpiReport<<"> To "<<proc<<"\t"<<pTest<<"\n"; app->mpiReport.flush();
				MPIUNLOCK;
			}
		}
	}
	
	while (activeProcesses > 0) {
		sleep(3);		
		for (int proc=0; proc<procCount; proc++ ){
			if (proc != app->nodeID) {
				int idSize = 2048;
				char modelIDs[idSize];
				float pResult = 0;
				MPILOCK;
				if (MPI::COMM_WORLD.Iprobe(proc, MPI_RESULT)) {
					MPI::COMM_WORLD.Recv(&pTest, 1, MPI::INT, proc, MPI_RESULT);
					MPI::COMM_WORLD.Recv(&pResult, 1, MPI::FLOAT, proc, MPI_RESULT);
					MPI::COMM_WORLD.Recv(modelIDs, idSize, MPI::CHAR, proc, MPI_RESULT);

					app->AppendToMORDist(dist, pTest, modelIDs, pResult);


					pTest = app->GetNextPTest();
					if (pTest >= 0)  {
						MPI::COMM_WORLD.Send(&pTest, 1, MPI::INT, proc, MPI_TESTID);
						app->mpiReport<<"> To "<<proc<<"\t"<<pTest<<"\n"; app->mpiReport.flush();
					}
					else
						activeProcesses--;
				}
				MPIUNLOCK;
			}
		}
	}
	pTest = -1;
	MPILOCK;
	for (int proc=0; proc<procCount; proc++) {
		if (proc!=app->nodeID) {
			MPI::COMM_WORLD.Send(&pTest, 1, MPI::INT, proc, MPI_TESTID);
		}
	}
	MPIUNLOCK;
	return NULL;
}

void AppMdrPDT::RunPTests(Distribution::PTestDistribution* dist) {
	pthread_mutex_lock(&ptestLock);
	Evaluation::TStatistic eval(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());
	pthread_mutex_unlock(&ptestLock);

	int testNumber = GetNextPTest();

	while (testNumber > 0) {
		Utility::Random rnd(seed+testNumber);
		GenotypeRepository test;

		//Create the test data partition
		pthread_mutex_lock(&ptestLock);
		repo.InitializeData(test, rnd, xvCount, true);
		pthread_mutex_unlock(&ptestLock);

		vector<string> exclusions;
		if (cfg.appConfiguration->GetLines("EXCLUDE_LOCUS", exclusions))
			test.InitExclusionList(exclusions);
		Evaluation::TStatistic testEval(xvCount, repo.GetPedigreeCount(), test.GetIndividualCount());
		Evaluation::TFinalReport ptestWinner(Evaluation::EvaluationMethod::maxModelSize, xvCount, 1);

		testEval.BasicEval(&test, ptestWinner);
		//We use an omnibus distribution, so we need to get the best of all N orders
		string bestID = "";
		float bestScore = 0.0;
		for (int o=Evaluation::EvaluationMethod::minModelSize-1; o<Evaluation::EvaluationMethod::maxModelSize;o++) {
			string modelID;
			float curScore = ptestWinner.GetBestMOR(o, modelID);
			if (curScore>bestScore) {
				bestScore=curScore;
				bestID = modelID;
			}
		}
		AppendToMORDist(dist, testNumber, bestID.c_str(),bestScore);
		testNumber = GetNextPTest();
	}
}

void AppMdrPDT::AppendToMORDist(Distribution::PTestDistribution* dist, int testID, const char *model, float score) {
	if (nodeID != masterID) {
		int idSize = 2048;
		mpiReport<<">"<<testID<<"\t"<<model<<"\t"<<score<<"\n"; mpiReport.flush();
		MPI::COMM_WORLD.Send(&testID, 1, MPI::INT, masterID, MPI_RESULT);
		MPI::COMM_WORLD.Send(&score, 1, MPI::FLOAT, masterID, MPI_RESULT);
		MPI::COMM_WORLD.Send(model, idSize, MPI::CHAR, masterID, MPI_RESULT);
	}
	else { 
		pthread_mutex_lock(&ptestLock);
		dist->AppendTest(1, testID-1, model, score);
		pthread_mutex_unlock(&ptestLock);
		std::cout<<".";std::cout.flush();
		mpiReport<<"<"<<testID<<"\t"<<model<<"\t"<<score<<"\n"; mpiReport.flush();
		if (ptestShortCircuit > 0) {
			if (score>winningScore) 
				ptestExceptions++;
			if (ptestExceptions>ptestShortCircuit) {
				cerr<<"\n\nToo many permutations have outperformed the top model (Model: "<<
					winningID<<" MOR="<<winningScore<<".) Halting execution prematurely\n";
				ptestsRequested=ptestCount;
			}
		}
	}
}

void AppMdrPDT::MpiPrep()	{
	stringstream nodeLogFilename;
	MPILOCK;
	nodeID = MPI::COMM_WORLD.Get_rank();
	processCount = MPI::COMM_WORLD.Get_size();
	masterID = 0;
	MPIUNLOCK;
	nodeLogFilename<<cfg.GetFilename()<<".mpi."<<nodeID<<".log";
	mpiReport.open(nodeLogFilename.str().c_str(), ios::out);
	
	if (nodeID != masterID) {
		origCout = std::cout.rdbuf();
		std::cout.rdbuf(mpiReport.rdbuf());
	}
}

void AppMdrPDT::MpiClose() {
	MPI::Finalize();
	mpiReport.close();

	if (nodeID != masterID) {
		std::cout.rdbuf(origCout);
	}
}
#endif ///USE_MPI


}
