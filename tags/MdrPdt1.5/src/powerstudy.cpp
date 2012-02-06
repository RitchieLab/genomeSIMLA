//
// C++ Implementation: powerstudy
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
pthread_mutex_t output_lock = PTHREAD_MUTEX_INITIALIZER;
#define TESTLOCK pthread_mutex_lock(&input_lock)
#define TESTUNLOCK pthread_mutex_unlock(&input_lock)
#define REPORTLOCK pthread_mutex_lock(&output_lock)
#define REPORTUNLOCK pthread_mutex_lock(&output_lock)

#else
//We don't want these to mean anything if we arne't using MPI
#define TESTLOCK
#define TESTUNLOCK
#define REPORTLOCK
#define REPORTUNLOCK

#endif 	//USE_MPI



#include "powerstudy.h"
#include "powerrepos.h"
#include <boost/timer.hpp>
#include "evalbalancedaccuracypdt.h"
#ifndef TEST_APP
int main(int argc, char **argv) {
	MDR::Power::PowerStudy app;

	if (app.ParseCmdLine(argc, argv))
		app.Start();
	return 0;
}
	
#endif

namespace MDR {

namespace Power {

using namespace boost;

PowerStudy::PowerStudy()
{
	appname="mdr-pdt-power";
	appfunction="MDR-PDT Power Study";
	authors="Marylyn Ritchie & Eric Torstenson\nPlease forward any comments or errors to: mdr-pdt@chgr.mc.vanderbilt.edu";
	major=1;
	minor=0;
	configuration = EseConfiguration::Instance();
}


PowerStudy::~PowerStudy()
{
	EseConfiguration::Release();
}


bool PowerStudy::ParseCmdLine(int argc, char **argv) {
	bool success=false;
	if (argc < 4) 
		PrintHelp();
	else {
		configFile = argv[1];
	
		LineParser lp;
		lp.Parse(configFile.c_str(), configuration);
		success=configuration->Validate();
		
		if (!success)
			cout<<configuration->errorMsg;
		else 
			configuration->PostLoad();
		
		targetModel = argv[2];

		for (int i =3; i < argc; i++) {
			inputfiles.push_back(argv[i]);
		}
	}
	return success;
}


bool PowerStudy::VerifyConfiguration() {
	bool isGood = configuration->analysisStyle == EseConfiguration::PDT;

	return isGood;
}
/**
 * Prints the help contents
 */
void PowerStudy::PrintHelp()	{
	PrintBanner();
	cout<<"usage: "<<appname<<" <configuration file> model_id data_file [data_file ...]\n";
	cout<<"\nThis application is intended to allow the user to calculate the power of a set of configuration parameters\n";
	cout<<"* A power study is performed on a single configuration file with multiple data files. \n";
	cout<<"* The results are tabulated and reported at the end. \n";
	cout<<"* Certain details (like p-tests) are ignored in the configuration file at this time. \n";
}

bool PowerStudy::NextDataset(string &filename) {
	TESTLOCK;
	bool doContinue = fileIdx<inputfiles.size();

	if (doContinue)
		filename = inputfiles[fileIdx++];

	TESTUNLOCK;

	return doContinue;
}
/**	
 * Starts execution
 */
void PowerStudy::Start()	{
	PrintBanner();

	timer progress;
	double totalTime = 0.0;

	configuration->SetFilename( inputfiles[0].c_str());
	configuration->GenerateReport(cout);

	string filename;

#ifdef USE_MPI
	pthread_t watcherThread; 
	threadArgs args(this, dist);	
	pthread_create(&watcherThread, NULL, &ProcessWatcher, (void*)&args);
#endif

	while (NextDataset(filename)) {
		totalTime += PedigreeSearch(filename.c_str());
	}

	cout<<"\n\n\n";
	cout<<setw(30)<<right<<"Model ID       "<<"     Frequency\n";
	
	map<string, uint>::iterator itr = power.begin();
	map<string, uint>::iterator end   = power.end();
	for (; itr != end; itr++) 
		cout<<setw(30)<<right<<itr->first<<setw(8)<<itr->second<<"\n";

	cout<<"\n\n--------------------------------\n";
	cout<<setw(25)<<right<<"Target         "<<setw(10)<<" "<<setw(10)<<"Average"<<setw(10)<<"Average"<<setw(10)<<"Time   "<<endl;
	cout<<setw(25)<<right<<"Model            "<<setw(10)<<"Power"<<setw(10)<<"T   "<<setw(10)<<"Time  "<<setw(10)<<"Expired"<<endl;
	cout<<setw(25)<<right<<targetModel;
	cout<<setw(10)<<right<<setprecision(2)<<((float)((power[targetModel].freq)/(float)datafileCount) * 100.0);
	cout<<setw(10)<<right<<setprecision(2)<<power[targetModel].GetAverageT();
	cout<<setw(10)<<right<<totalTime/(float)datafileCount;
	cout<<setw(10)<<right<<progress.elapsed()<<"\n";
	cout<<"--------------------------------\n";
		
}


double PowerStudy::PedigreeSearch(const char *filename) {
	timer progress;
	configuration->SetFilename(filename);
	SnpRepository 			snps;				///<The original set of snps

	uint foldCount = configuration->crossValidationCount;

	//Acquire the input parser from the configuration
	GtFileParser *fileParser=configuration->GetInputFileParser(true);
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );

	if (snps.GetSnpCount() == 0) {
		cout<<"\nNo data to search over. Please check that the file that you specified exists\n";
		abort();
	}
	

	//Close the locus log, if it exists
	if (configuration->logLocus)
		configuration->logLocus->Close();

	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	//Create a report to store the top models
	PowerRepos bestScore(configuration->reportModelCount);
	SnpEvalSuite *e = configuration->GetPdtEval( false, 0 );
	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)e->GetTrainer();

	//The distribution object needs to be passed to the various ptests
	PermutationTestDist *dist = e->GetDistribution();
	bestScore.distribution = dist;
	BasicLog *reportLog = configuration->GetReportLog();
	double evalTime;
	//If some analysis is required, do it now
	if (e) {
		double loadTime = progress.elapsed();
		
		progress.restart();
		//Perform the search

		uint modelCount = snps.Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, e);
		//cout<<"\n\nPedigree Search of "<<e->GetSnpCount()<<" : "<<modelCount<<" completed in "<<progress.elapsed()<<" seconds\n";
		
		evalTime = progress.elapsed();
	
		bestScore.Sort();
		PowerRepos::ReportModel *model = bestScore.GetReportModel(0, 0);
		string topModelID = model->GetLabel();
		SnpAligned *curModel=snps.GetSnp(topModelID.c_str());
		power[topModelID].freq++;
		power[topModelID].totalT += curModel->GetLastMdEval();

		uint reportCount = 5;
		if (bestScore.GetEntryCount(0) < reportCount)
			reportCount = bestScore.GetEntryCount(0);
	
		for (uint i=0; i<reportCount; i++) {
			model = bestScore.GetReportModel(0, i);			//We don't care about loci, since our repository does't descriminate
			curModel = snps.GetSnp(model->GetLabel().c_str());

			for (uint fold=0; fold<foldCount; fold++) {
				stringstream details;
				cout<<setw(5)<<" "<<fileIdx<<setw(25)<<filename<<setw(5)<<i<<setw(8)<<setprecision(4)<<loadTime<<setw(8)<<setprecision(4)<<evalTime;

				evalPdt->EvaluateVerbose(fold, curModel, details);
			}
		}
		
		delete e;
		
	}
	return evalTime;
}



#ifdef USE_MPI 

struct threadArgs {
	ModelFinder *main;
	PermutationTestDist *fDist;
	threadArgs(ModelFinder *main, PermutationTestDist *fDist) : main(main), fDist(fDist) {}
};

void *PowerStudy::ProcessWatcher(void *arg) {
	ModelFinder *main = ((threadArgs*)arg)->main;
	PermutationTestDist *dist = ((threadArgs*)arg)->fDist;

	int procCount = main->processCount;
	int maxModelSize = main->configuration->comboEnd;
	int activeProcesses = 0;
	int tests[2];	
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {


			if (NextDataset(filename)) {
				activeProcesses++;
				MPI::COMM_WORLD.Isend(filename.c_str(), strlen(filename.c_str(), MPI::CHAR, proc, TESTID);
				
			}
		}
	}
	cout<<"#\tMaster Node, "<<main->nodeID<<", initiated "<<activeProcesses<<" p-tests\n";

	//Wait for the recipients to finish and respond via round robin/
	//This thread will NEVER stop if a node can't return after being assigned
	//a legitimate test ID
	while (activeProcesses > 0) {
		/*  I need to figure out how to get what I need from the power results
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
					
					cout<<"#$#$#$ "<<idSize<<" - "<<modelIDs<<"\n";
					
				
					main->AppendPTest(tests[0], tests[1], pResults, dist, modelIDs);
					if (main->NextTestID(tests[0])) {
						MPI::COMM_WORLD.Isend(tests, 1, MPI::INT, proc, TESTID);
					}
					else
						activeProcesses--;
				}
			}
		}
		*/
	}

	//Signal each thread to go away
	tests[0] = -1;
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			cout<<"#\tAttempting to shut down node "<<proc<<"\n";
			MPI::COMM_WORLD.Send(tests, 1, MPI::INT, proc, TESTID);
		}
	}
	return NULL;
	
}
#endif

}

}
