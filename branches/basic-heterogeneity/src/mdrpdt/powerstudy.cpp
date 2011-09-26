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




#include "powerstudy.h"
#include <boost/timer.hpp>
#include "evalbalancedaccuracypdt.h"
#include "genetics/snprepostxtsorted.h"



int main(int argc, char **argv) {
	MDR::Power::PowerStudy app;

	if (app.ParseCmdLine(argc, argv) && app.PerformSanityCheck())
		app.Start();
	return 0;
}
	

namespace MDR {

namespace Power {


#define RESULTS_FILENAME 			0
#define RESULTS_XV_CONSISTENCY 		1
#define RESULTS_LABEL 				2
#define RESULTS_STATS 				3
#define RESULTS_META 				4
#define RESULTS_TESTID				5
#define LOAD_DISTRIBUTIONS			6
#define DISTRIBUTIONS_LOADED		7

using namespace boost;

PTestDistribution *PowerStudy::ModelHolder::fDist = NULL;
PTestDistribution *PowerStudy::ModelHolder::oDist = NULL;
PTestDistribution *PowerStudy::ModelHolder::pDist = NULL;
bool PowerStudy::VerbosePower = false;


PowerStudy::PowerStudy() : fileIdx(0), useReferenceDistribution(true)
{
#ifdef USE_MPI
	appname="parallel mdr-pdt-power";
#else
	appname="mdr-pdt-power";
#endif
	appfunction="MDR-PDT Power Study";
	authors="Marylyn Ritchie & Eric Torstenson\nPlease forward any comments or errors to: mdr-pdt@chgr.mc.vanderbilt.edu";
	major=1;
	minor=2;
	bugFixes=5;
	configuration = EseConfiguration::Instance();
	rand = &Utility::Random::globalGenerator;

}


PowerStudy::~PowerStudy()
{
	EseConfiguration::Release();
}


bool PowerStudy::ParseCmdLine(int argc, char **argv) {
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
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

	//If we are the master, we should grab the file from the list. Otherwise, it's coming over the network
	if (nodeID == masterID) {
	
		if (doContinue)
			filename = inputfiles[fileIdx++];
	}
#ifdef USE_MPI
	else {
		char fn[4096];
		int n=4096;
	
		MPI::COMM_WORLD.Recv(fn, n, MPI::CHAR, masterID, RESULTS_FILENAME);
		doContinue = !(fn==NULL || strlen(fn)==0);
		if (doContinue)
			filename=fn;
	}			
#endif	
	TESTUNLOCK;

	return doContinue;
}

#ifdef USE_MPI 

struct threadArgs {
	PowerStudy *main;
	Distribution::PTestDistribution *fDist;
	threadArgs(PowerStudy *main, Distribution::PTestDistribution *fDist) : main(main), fDist(fDist) {}
};

void *PowerStudy::ProcessWatcher(void *arg) {
	PowerStudy *main = ((threadArgs*)arg)->main;
	//Distribution::PTestDistribution *dist = ((threadArgs*)arg)->fDist;

	int procCount = main->processCount;
	int maxModelSize = main->configuration->comboEnd;
	int activeProcesses = 0;
	int tests[2];	
	string filename;
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			if (main->NextDataset(filename)) {
				activeProcesses++;
				MPI::COMM_WORLD.Isend(filename.c_str(), strlen(filename.c_str()) + 1, MPI::CHAR, proc, RESULTS_FILENAME);
				main->log<<proc<<" -> "<<filename<<"\n";
				cout<<"Master Node Sent "<<proc<<" the file, "<<filename<<"\n";
			}
		}
	}
	main->log<<"Master Node, "<<main->nodeID<<", initiated "<<activeProcesses<<" p-tests\n";

	//Wait for the recipients to finish and respond via round robin/
	//This thread will NEVER stop if a node can't return after being assigned
	//a legitimate test ID
	while (activeProcesses > 0) {

		for (int proc=0; proc<procCount; proc++) {
			if (proc!=main->nodeID) {		
				if (MPI::COMM_WORLD.Iprobe(proc, RESULTS_FILENAME)) {
					char fn[4096];
					for (int i=0; i<maxModelSize; i++) {
						char xvCons[24];
						char label[128];
						float stats[5];
						int denom[5];
						int meta[4];

						MPI::COMM_WORLD.Recv(fn, 4096, MPI::CHAR, proc, RESULTS_FILENAME);
						MPI::COMM_WORLD.Recv(xvCons, 24, MPI::CHAR, proc, RESULTS_XV_CONSISTENCY);
						MPI::COMM_WORLD.Recv(label, 128, MPI::CHAR, proc, RESULTS_LABEL);
						MPI::COMM_WORLD.Recv(stats, 5, MPI::FLOAT, proc, RESULTS_STATS);
						MPI::COMM_WORLD.Recv(denom, 5, MPI::INT, proc, RESULTS_STATS);
						MPI::COMM_WORLD.Recv(meta, 4, MPI::INT, proc, RESULTS_META);

														
						assert(meta[1] == i);

						ModelStatistics st;
						st.sumTrainingStatistic=stats[0];
						st.sumClassError=stats[1];
						st.sumTestingStatistic=stats[2];
						st.sumPredError=stats[3];
						st.sumMOR=stats[4];
						st.xvConsistency=xvCons;
						st.xvConst = meta[3];
						st.label=label;

						st.trainingCount=denom[0];
						st.clsCount=denom[1];
						st.testingCount=denom[2];
						st.predCount=denom[3];
						st.orCount=denom[4];

						//testID = meta[0];
						st.foldCount = meta[2];
						//modelSize = i;


						
						main->SaveResults(proc, fn, meta[0], meta[1], st);
					}				

					if (main->NextDataset(filename)) {	
						strcpy(fn, filename.c_str());
						main->log<<proc<<" -> "<<fn<<"\n";
						MPI::COMM_WORLD.Isend(fn, 4096, MPI::CHAR, proc, RESULTS_FILENAME);
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
			main->log<<"Attempting to shut down node "<<proc<<"\n";
			char fn[4096];
			fn[0]='\0';
			MPI::COMM_WORLD.Isend(fn,4096, MPI::CHAR, proc, RESULTS_FILENAME);
		}
	}

	return NULL;
	
}


void *PowerStudy::PTestWatcher(void *arg) {

	PowerStudy *main = ((threadArgs*)arg)->main;
	//Distribution::PTestDistribution *dist = ((threadArgs*)arg)->fDist;
	//We only want to do this if we are building a reference distribution
	if (!main->useReferenceDistribution)
		return NULL;

	int procCount = main->processCount;
	int maxModelSize = main->configuration->comboEnd;
	int activeProcesses = 0;
	int tests[2];	
	size_t testID;
	//Send off the first batch of tests via non-blocking sends
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			if (main->NextTestID(testID)) {
				activeProcesses++;

				/**
	    		 * EST - This is really ugly, but I need to figure out why testID is size_t
	 			 */
				uint curID = (uint)testID;
				MPI::COMM_WORLD.Send(&curID, 1, MPI::INT, proc, TESTID);
				main->log<<proc<<" -> "<<testID<<"\n";

			}
		}
	}
	main->log<<"Master Node, "<<main->nodeID<<", initiated "<<activeProcesses<<" p-tests\n";

	//Wait for the recipients to finish and respond via round robin/
	//This thread will NEVER stop if a node can't return after being assigned
	//a legitimate test ID
	while (activeProcesses > 0) {
		for (int proc=0; proc<procCount; proc++) {

			vector<ModelStatistics> recStatistics;

			if (proc!=main->nodeID) {		
				if (MPI::COMM_WORLD.Iprobe(proc, RESULTS_TESTID)) {
					char fn[4096];
					float predictedT=0.0;
					uint modelSize = 0;
					for (int i=0; i<maxModelSize; i++) {
						ModelStatistics st;

						char xvCons[24];
						char label[128];
						float stats[5];
						int denom[5];
						int meta[4];

						MPI::COMM_WORLD.Recv(&testID, 1, MPI_INT, proc, RESULTS_TESTID);
						//MPI::COMM_WORLD.Recv(fn, 4096, MPI::CHAR, proc, RESULTS_FILENAME);
						MPI::COMM_WORLD.Recv(xvCons, 24, MPI::CHAR, proc, RESULTS_XV_CONSISTENCY);
						MPI::COMM_WORLD.Recv(label, 128, MPI::CHAR, proc, RESULTS_LABEL);
						MPI::COMM_WORLD.Recv(stats, 5, MPI::FLOAT, proc, RESULTS_STATS);
						MPI::COMM_WORLD.Recv(denom, 5, MPI::INT, proc, RESULTS_STATS);
						MPI::COMM_WORLD.Recv(meta, 4, MPI::INT, proc, RESULTS_META);

	
						assert(meta[1] == i);
	
					
						st.sumTrainingStatistic=stats[0];
						st.sumClassError=stats[1];
						st.sumTestingStatistic=stats[2];
						st.sumPredError=stats[3];
						st.sumMOR=stats[4];
						st.xvConsistency=xvCons;
						st.xvConst = meta[3];
						st.label=label;

						st.trainingCount=denom[0];
						st.clsCount=denom[1];
						st.testingCount=denom[2];
						st.predCount=denom[3];
						st.orCount=denom[4];

						testID = meta[0];
						st.foldCount = meta[2];
						modelSize = i;

						recStatistics.push_back(st);

						//cout<<"Node "<<proc<<" Test: "<<testID<<" -> "<<st.label<<" "<<st.xvConst<<" "<<st.sumTestingStatistic<<" "<<st.sumPredError<<" "<<st.sumMOR<<"\n";

					}				
//					cout<<"Saving Results for test#"<<testID;
					main->SavePTests(proc, testID, recStatistics);
//					cout<<"....Saved!\n";
//					if (predictedT > 0.0)
//						main->SavePTest(proc, testID, modelSize, st);


					if (main->NextTestID(testID)) {
						cout<<"Sending new assignment ("<<testID<<") to node "<<proc;
						MPI::COMM_WORLD.Send(&testID, 1, MPI::INT, proc, TESTID);
						main->log<<proc<<" -> "<<testID<<"\n";


						cout<<"....Done!\n";

					}
					else{
						cout<<"We are done- gotta wait for the rest of the nodes to finish up!\n";
						activeProcesses--;
					}
				}
			}
		}
	}

	cout<<"OK, sending the children message to stop waiting for PTests\n";
	//Signal each thread to go away
	for (int proc=0; proc<procCount; proc++) {
		//Make sure we aren't telling ourselves to do anything
		if (proc != main->nodeID) {
			main->log<<"Attempting to shut down node "<<proc<<"\n";
			int msg=0;
			MPI::COMM_WORLD.Isend(&msg,1, MPI::INT, proc, TESTID);
		}
	}

	cout<<"Message sent to the "<<procCount-1<<" nodes\n";
	return NULL;
	
}
#endif


/**	
 * Starts execution
 */
void PowerStudy::Start()	{
	REPORTLOCK;							//We are printing stuff in the recording, so let's make sure that doesn't happen here

	MpiPrep();
	
	stringstream nodeLogFilename;
	nodeLogFilename<<configuration->reportName;
	nodeLogFilename<<".node"<<nodeID<<".log";
	log.open(nodeLogFilename.str().c_str(),ios::out|ios::trunc);
	log<<"Node "<<nodeID<<" up!\n";

	if (nodeID == masterID)
		PrintBanner();

	timer progress;
	double totalTime = 0.0;

	configuration->SetFilename( inputfiles[0].c_str());

	if (nodeID == masterID) 
		configuration->GenerateReport(cout);

	string filename;


	

#ifdef USE_MPI
	pthread_t watcherThread; 
	threadArgs args(this, NULL);
	//This is only appropriate if we are the master node
	if (nodeID == masterID) {
		cout<<setw(45)<<right<<"Total MPI Node Count: "<<processCount<<endl;

		pthread_create(&watcherThread, NULL, PTestWatcher, (void*)&args);
		log<<"Started PTest watcher thread\n";

	}
#endif
	REPORTUNLOCK;
	//Do the PTests 
	size_t testID = 0;
	int modelCount = 0;

	rand->Seed(configuration->pTestSeed);

	SnpRepository 			snps;				///<The original set of snps	

	//Acquire the input parser from the configuration
	fileParser=configuration->GetInputFileParser(true);

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );

	configuration->ValidateSnpCount(snps.GetSnpCount());

	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	//Create a report to store the top models
	SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);

	SnpEvalSuite *e = configuration->GetPdtEval( false, 0 );
//	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)e->GetTrainer();

	//The distribution object needs to be passed to the various ptests
	fDist = e->GetFitnessDist();
	ModelHolder::fDist = fDist;

	oDist = e->GetOddsRatioDist();
	ModelHolder::oDist = oDist;

	pDist = e->GetPEDist();
	ModelHolder::pDist = pDist;

	bestScore.distribution = fDist;

	//Build up the reference distribution
	if (useReferenceDistribution) {
		if (nodeID == masterID) {
			cout<<"Setting up referential distributions ("<<configuration->pTestCount<<")\n";	
			cout<<"Referential distributions based on dataset: "<< inputfiles[0].c_str()<<"\n";
		}
		while (NextTestID(testID)) {
			vector<ModelStatistics> bestModels;
			RunPTest(testID, modelCount, bestModels);
			SavePTests(nodeID, testID - 1, bestModels);
	
		}
#ifdef USE_MPI
		if (nodeID == masterID) {
			cout<<"Waiting for the PTest Watcher thread to join";
			log<<"Waiting for the watcher thread to join\n";
			pthread_join(watcherThread, NULL);
			cout<<"....completed\n"; 
			log<<"Completed.\n";
		}
#endif
		//If we are the Master node, we have the distribution available to save to file.
		//The other nodes need to load it. 
		ManageReferenceDist();
	}

	cout<<"Getting ready to Start the regular search!\n";
	REPORTLOCK;
#ifdef USE_MPI
	if (nodeID == masterID) {
		cout<<"Initiating the watcher thread for results\n";
		pthread_create(&watcherThread, NULL, &ProcessWatcher, (void*)&args);
		log<<"Master watcher thread created and running\n";
	}
#endif
	REPORTUNLOCK;
	while (NextDataset(filename)) {
		totalTime += PedigreeSearch(filename.c_str(), 0, 0);
	}

	if (nodeID == masterID) {
#ifdef USE_MPI
		log<<"Waiting for the watcher thread to join\n";
		pthread_join(watcherThread, NULL); 
		log<<"Completed. We can now report the final results\n";
		//((Genetics::Distribution::PTestDistribution*)fDist)->Report(&cout);
		//fDist->Report(&cout);
		//oDist->Report(&cout);
#endif

		cout<<"\n\n\n";

		cout<<"Power based on T-Statistic: \n";
		ReportPower(power);

		cout<<"\n\nPower Based on Prediction Error\n";
		ReportPower(pePower);

		cout<<"\n\nPower Based on Matched Odds Ration\n";
		ReportPower(orPower);

	}
	MpiClose();
}

void PowerStudy::ReportPower(map<string, ModelResults>& power) {
		cout<<setw(30)<<right<<"Model ID       "<<setw(10)<<"Frequency"<<setw(10)<<"Avg T"<<setw(10)<<"Avg P.E."<<setw(10)<<"Avg MOR"<<"\n";
		
		map<string, ModelResults>::iterator itr = power.begin();
		map<string, ModelResults>::iterator end   = power.end();
		for (; itr != end; itr++) 
			cout<<setw(30)<<right<<itr->first<<setw(10)<<itr->second.freq
				<<setw(10)<<itr->second.GetAverageTesting()
				<<setw(10)<<itr->second.GetAveragePrd()
				<<setw(10)<<itr->second.GetAverageOR()<<"\n";
	

		int datafileCount = inputfiles.size();
		cout<<"\n\n-----------------------------------------------------------------";
		if (configuration->crossValidationCount > 1)
			cout<<"----------------------------";
		cout<<endl;

		cout<<setw(35)<<""<<setw(10);
		cout<<"Avg  "<<setw(10)<<"Avg  ";
		if (configuration->crossValidationCount > 1)
			cout<<setw(10)<<"Avg  "<<setw(10)<<"Avg  ";
		cout<<setw(10)<<"Matched"<<endl;


		cout<<setw(25)<<right<<"Target";
		cout<<setw(10)<<" "<<setw(10)<<"Class."<<setw(10)<<"Train.";
		if (configuration->crossValidationCount > 1)
			cout<<setw(10)<<"Pred."<<setw(10)<<"Testing";
		cout<<setw(10)<<"Odds";
		cout<<endl;

		cout<<setw(25)<<right<<"Model"<<setw(10)<<"Power";
		cout<<setw(10)<<"Error"<<setw(10)<<"MDR-PDT";
		if (configuration->crossValidationCount > 1)
			cout<<setw(10)<<"Error"<<setw(10)<<"MDR-PDT"<<setw(10);
		cout<<setw(10)<<"Ratio"<<endl;
		//cout<<"Time  "<<setw(10);
		//cout<<"Expired"<<endl;
		cout<<setw(25)<<right<<targetModel;
		float thePower = (float)power[targetModel].freq/(float)datafileCount * 100.0;
		cout<<setw(10)<<setprecision(2)<<thePower;

		
		if (configuration->crossValidationCount == 1 ) {
			cout<<setw(9)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<power[targetModel].GetAveragePrd()<<"%";
			cout<<setw(10)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<power[targetModel].GetAverageTraining();
		} else {
			cout<<setw(9)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<power[targetModel].GetAverageCls()<<"%";
			cout<<setw(10)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<power[targetModel].GetAverageTraining();
			cout<<setw(9)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<power[targetModel].GetAveragePrd()<<"%";
			cout<<setw(10)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<power[targetModel].GetAverageTesting();
		}
		cout<<setw(10)<<right<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<power[targetModel].GetAverageOR();
		//cout<<setw(10)<<right<<totalTime/(float)datafileCount;
		//cout<<setw(10)<<right<<progress.elapsed();
		cout<<"\n-----------------------------------------------------------------";
		if (configuration->crossValidationCount > 1)
			cout<<"----------------------------";
		cout<<endl;

}

void PowerStudy::ManageReferenceDist() {
	if (nodeID == masterID) 
	{

		cout<<"Managing the reference distributions\n";
		log<<"Referential Distributions: \n\n";
	
		log<<"Fitness Distribution\n";
		fDist->DumpDistribution(&cout);
	
		log<<"\n\nMOR Distribution\n";
		oDist->DumpDistribution(&cout);
	
		log<<"\n\nPE Distribution\n";
		pDist->DumpDistribution(&cout);
	}

#ifdef USE_MPI
	string morRefDist = configFile + string(".mordist");
	string peRefDist  = configFile + string(".pedist");
	string fRefDist   = configFile + string(".fdist");
	int data;				///<Boring buffer used by MPI

	if (nodeID == masterID) 
	{

		cout<<"Managing the reference distributions\n";
		log<<"Referential Distributions: \n\n";
	
		log<<"Fitness Distribution\n";
		fDist->DumpDistribution(&cout);
	
		log<<"\n\nMOR Distribution\n";
		oDist->DumpDistribution(&cout);
	
		log<<"\n\nPE Distribution\n";
		pDist->DumpDistribution(&cout);
		
		log.flush();

		//Save the distributions to files
		ofstream distfile(morRefDist.c_str(), ios_base::out);
		oDist->DumpDistribution(&distfile);
		distfile.close();
		
		distfile.open(peRefDist.c_str(), ios_base::out);
		pDist->DumpDistribution(&distfile);
		distfile.close();
		
		distfile.open(fRefDist.c_str(), ios_base::out);
		fDist->DumpDistribution(&distfile);
		distfile.close();

		//Iterate over each of the slaves and tell when it's time to load distributions
		for (int proc=0; proc<processCount; proc++) {
			//Make sure we aren't telling ourselves to do anything
			if (proc != nodeID) {
				cout<<"Signaling node "<<proc<<" to load distributions\n";
				MPI::COMM_WORLD.Isend(&data, 1, MPI::INT, proc, LOAD_DISTRIBUTIONS);
				cout<<"Waiting for acknowledgement of distributions loaded\n";
				MPI::COMM_WORLD.Recv(&data, 1, MPI::INT, proc, DISTRIBUTIONS_LOADED);
				log<<proc<<" Instructed to load distribution\n";
			}
		}
		log<<"Master Node, "<<nodeID<<" issued load distribution to all nodes\n";
		cout<<"Master Node "<<nodeID<<" issued load distribution to all nodes\n";
	}
	else 
	{

		//Wait for the message- otherwise, we might try to load an invalid set
		cout<<"Node "<<nodeID<<" awaiting cue to load distributions\n";
	    MPI::COMM_WORLD.Recv(&data, 1, MPI::INT, masterID, LOAD_DISTRIBUTIONS);
		cout<<"Got it!\n";

		cout<<"Node "<<nodeID<<" trying to load file: "<<morRefDist<<"\n";
		ifstream distfile(morRefDist.c_str(), ios_base::in);
		cout<<"File open\n";
		oDist->LoadDistribution(&distfile);
		cout<<"File loaded\n";
		distfile.close();
		cout<<"file closed\n";

		cout<<"Node "<<nodeID<<" trying to load file: "<<peRefDist<<"\n";
		distfile.open(peRefDist.c_str(), ios_base::in);
		pDist->LoadDistribution(&distfile);
		distfile.close();
		cout<<"file closed\n";

		cout<<"Node "<<nodeID<<" trying to load file: "<<fRefDist<<"\n";
		distfile.open(fRefDist.c_str(), ios_base::in);
		fDist->LoadDistribution(&distfile);
		distfile.close();

		cout<<"All done. Acknowledgement sent\n";
		//Send acknowledgement
		MPI::COMM_WORLD.Send(&data, 1, MPI::INT, masterID, DISTRIBUTIONS_LOADED);

	}
	
#endif

}


/**
 * @brief Perform the actual search
 * @param filename - This is the name of the datafile to be used
 * @param seed 	   - This is the seed to be used for this set
 * @param testID   - This is a number from 1 to N if it's a p-test, otherwise, it should be a 0
 */
double PowerStudy::PedigreeSearch(const char *filename, uint seed, uint testID) {

	timer progress;
	configuration->SetFilename(filename);
	SnpRepository 			snps;				///<The original set of snps

	//uint foldCount = configuration->crossValidationCount;
	rand->Seed(configuration->pTestSeed + seed);

	//Acquire the input parser from the configuration
	fileParser=configuration->GetInputFileParser(true);

	//Load the data into the buffered parser and write to the locus log any snps that are discarded (and why)
	snps.ParseInputFile( fileParser, configuration->logLocus );
	configuration->ValidateSnpCount(snps.GetSnpCount());
	
	if (snps.GetSnpCount() == 0) {
		cout<<"\n"<<nodeID<<"("<<filename<<")  has no data to search over. Please check that the file that you specified exists\n";
		abort();
	}
	

	//Close the locus log, if it exists
	if (configuration->logLocus)
		configuration->logLocus->Close();

	//Get the overall status mask 
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	//Create a report to store the top models
	SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);

	if (!useReferenceDistribution)
		configuration->ResetDistributions();
	SnpEvalSuite *e = configuration->GetPdtEval(testID!=0, testID);
	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)e->GetTrainer();

	

	double evalTime;
	ModelStatistics topStatistics[configuration->comboEnd];
	//If some analysis is required, do it now
	if (e) {
		//double loadTime = progress.elapsed();
		
		progress.restart();
		//Perform the search
		snps.Evaluate(configuration->comboStart - 1, configuration->comboEnd - 1, &bestScore, e);
		//cout<<"\n\nPedigree Search of "<<e->GetSnpCount()<<" : "<<modelCount<<" completed in "<<progress.elapsed()<<" seconds\n";


		evalTime = progress.elapsed();

		//This just grabs the top model. 
		for (uint modelSize= configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
			SnpAligned* bestModel = NULL;
			topStatistics[modelSize] = e->GetTopModel(modelSize, bestModel);
			topStatistics[modelSize].GatherTestStatistics(&snps, evalPdt, NULL);


			SaveResults(nodeID, filename, 0, modelSize, topStatistics[modelSize]);
		}
	}

		if (VerbosePower)
			evalPdt->DisplayAnalyses(true);

	//Unless we are building individual distributions, we don't need to do this
	if (!useReferenceDistribution) {


		//Build up the distributions
		size_t i =0; 
		while (NextTestID(i)) {
			int modelCount;
			vector<ModelStatistics> topModels;
			//float searchTime = 
			RunPTest(i, modelCount, topModels);
	
			uint stIdx = 0;
	
			for (uint modelSize= 0; modelSize < configuration->comboEnd; modelSize++) {
				if (modelSize >= configuration->comboStart-1) {
					//If this fails the RunPTest probably failed
					assert(topModels.size() > stIdx);
					ModelStatistics &st = topModels[stIdx++];
					st.GatherTestStatistics(&snps, evalPdt, NULL);
					cout<<i<<"\tModel: "<<st.GetLabel()<<"\t";
					fDist->AppendTest(st.GetAvgTesting(), modelSize, st.GetLabel().c_str(), i);
					cout<<"\tFitness: "<<st.GetAvgTesting();
					oDist->AppendTest(st.GetOddsRatio(), modelSize, st.GetLabel().c_str(), i);
					cout<<"\tMOR: "<<st.GetOddsRatio();
					pDist->AppendTest(st.GetAvgPredictionError(), modelSize, st.GetLabel().c_str(), i);
					cout<<"\tPE: "<<st.GetAvgPredictionError()<<"\n";
				}
			}			
		}
	}

	//If some analysis is required, do it now
	if (e) {
		//double loadTime = progress.elapsed();
		
		string topORModel;
		string topPEModel;
		
		float bestOR = 1.1;
		float bestPE = 1.1;

		ModelStatistics topOR;
		ModelStatistics topPE;
		
		//This just grabs the top model. 
		for (uint modelSize= configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
			ModelStatistics &st = topStatistics[modelSize];

/*			cout<<"\n"<<modelSize<<" "<<setw(15)<<st.GetLabel()
				<<setw(10)<<st.GetAvgTesting()<<setw(10)<<st.GetAvgPredictionError()
				<<setw(10)<<st.GetOddsRatio()<<"\n";
	*/	

			float orPValue = 0.0;
			if (oDist)
				orPValue = oDist->GetPValue(st.GetOddsRatio(), modelSize);

			float pePValue = 0.0;
			if (pDist)
				pePValue = pDist->GetPValue(st.GetAvgPredictionError(), modelSize);

			if (orPValue < bestOR) {
				bestOR = orPValue;
				topORModel = st.GetLabel();
				topOR = st;
			}

			if (pePValue < bestPE) {
				bestPE = pePValue;
				topPEModel = st.GetLabel();
				topPE = st;
			}
				
		}

		ModelResults &r = pePower[topPEModel];
		r+=topPE;
		ModelResults &s = orPower[topORModel];
		s+=topOR;

		delete e;

	}
	return evalTime;
}

void PowerStudy::SavePTests(int originatorID, int testID, vector<ModelStatistics> bestModels) {
	uint modelCount = bestModels.size();

	//cout<<"Saving PTests "<<testID<<" From node: "<<originatorID<<"\n";

	//This will be used to discriminate between xv consistency & t-values 
	ModelStatistics topModel;
	for (uint i=0; i<modelCount; i++) {
		ModelStatistics &cur=bestModels[i];

		if (nodeID == masterID) {
//			cout<<"Writing to fitness distribution ("<<i<<")\n";
			if (fDist)
				fDist->AppendTest(cur.GetAvgTesting(), i, cur.GetLabel().c_str(), testID);
			if (topModel < cur)
				topModel = cur;

			//cout<<"PTest ("<<originatorID<<":"<<testID<<" o:"<<i<<") "<<" "<<cur.GetLabel()<<" "<<cur.GetAvgTesting()<<" "<<cur.GetAvgPredictionError()<<" "<<cur.GetOddsRatio()<<"\n";
		}
		else {
//			cout<<"Sending them to the Master Node\n";
			//This is the only place we should call this function
			SavePTest(originatorID, testID, i, cur);
		}
	}

	if (nodeID == masterID) {
//		cout<<"Writing to other distributions\n";

		//cout<<"Winning model selected: ("<<originatorID<<":"<<testID<<" o) "<<" "<<topModel.GetLabel()<<" "<<topModel.GetAvgTesting()<<" "<<topModel.GetAvgPredictionError()<<" "<<topModel.GetOddsRatio()<<"\n";
		cout<<".";
		cout.flush();
		if (oDist)
			oDist->AppendTest(topModel.GetOddsRatio(), 0, topModel.GetLabel().c_str(), testID);
		if (pDist)
			pDist->AppendTest(topModel.GetAvgPredictionError(), 0, topModel.GetLabel().c_str(), testID);
	}
}
void PowerStudy::SavePTest(int originatorID, int testID, int modelSize, ModelStatistics &st) {
	REPORTLOCK;

	if (nodeID == masterID) {
		assert(0);
		cout<<"\t#"<<originatorID<<"\tTest #"<<testID<<"\t"<<st.GetAvgTesting()<<":"
			<<st.GetAvgPredictionError()<<":"<<st.GetOddsRatio()<<"\n";
		if (oDist)
			oDist->AppendTest(st.GetOddsRatio(), modelSize, st.GetLabel().c_str(), testID);
		if (pDist)
			pDist->AppendTest(st.GetAvgTesting(), modelSize, st.GetLabel().c_str(), testID);
	}			
	else {
#ifdef USE_MPI
		float stats[5];
		int denom[5];
		char label[128];
		char consistency[24];
		char fn[4096];
		stats[0]=st.sumTrainingStatistic;
		stats[1]=st.sumClassError;
		stats[2]=st.sumTestingStatistic;
		stats[3]=st.sumPredError;
		stats[4]=st.sumMOR;

		denom[0]=st.trainingCount;
		denom[1]=st.clsCount;
		denom[2]=st.testingCount;
		denom[3]=st.predCount;
		denom[4]=st.orCount;

		int meta[4];
		meta[0]=testID;
		meta[1]=modelSize;
		meta[2]=st.foldCount;
		meta[3]=st.xvConst;

		//strcpy(fn,filename);
		

		strcpy(label,st.label.c_str());
		
	
		strcpy(consistency,st.xvConsistency.c_str());
		//log<<label<<" "<<consistency<<" "<<filename<<"\n";
		log<<label<<" "<<consistency<<"\n";


		MPI::COMM_WORLD.Send(&testID, 1, MPI::INT, masterID, RESULTS_TESTID);
		//MPI::COMM_WORLD.Recv(fn, 4096, MPI::CHAR, masterID, RESULTS_FILENAME);
		MPI::COMM_WORLD.Send(consistency, 24, MPI::CHAR, masterID, RESULTS_XV_CONSISTENCY);
		MPI::COMM_WORLD.Send(label, 128, MPI::CHAR, masterID, RESULTS_LABEL);
		MPI::COMM_WORLD.Send(stats, 5, MPI::FLOAT, masterID, RESULTS_STATS);
		MPI::COMM_WORLD.Send(denom, 5, MPI::INT, masterID, RESULTS_STATS);
		MPI::COMM_WORLD.Send(meta, 4, MPI::INT, masterID, RESULTS_META);
#endif
	}

	REPORTUNLOCK;
	
}

void PowerStudy::SaveResults(uint nID, const char *filename, uint position, uint modelSize, ModelStatistics &st) {
	static bool printedHeader = false;
	REPORTLOCK;

	if (nodeID == masterID) {
		ModelHolder v(filename, position, modelSize, st);
		if (!printedHeader) {
			v.GenerateHeader(&cout);
			printedHeader = true;
		}			


		string topModelID = st.label;
		log<<nID<<" <- ("<<filename<<") "<<topModelID<<" "<<st.GetAvgTraining()<<"\n";
		ModelResults &res=power[topModelID];
		res+=st;
		//power[topModelID].freq++;
		/*if (st.foldCount > 1)
			power[topModelID].totalT += st.GetAvgTraining();
		else
			power[topModelID].totalT += st.GetAvgTesting();
		*/
		v.GenerateReport(&cout);
		statistics.push_back(v);
	}
#ifdef USE_MPI 
	else	{
		//OK, we need to extract everything from the st and send it back to the Master node
		//For now, let's just send the details for the averages and other little tidbits
		float stats[5];	
		int denom[5];
		char label[128];
		char consistency[24];
		char fn[4096];
		stats[0]=st.sumTrainingStatistic;
		stats[1]=st.sumClassError;
		stats[2]=st.sumTestingStatistic;
		stats[3]=st.sumPredError;
		stats[4]=st.sumMOR;

		denom[0]=st.trainingCount;
		denom[1]=st.clsCount;
		denom[2]=st.testingCount;
		denom[3]=st.predCount;
		denom[4]=st.orCount;


		int meta[4];
		meta[0]=position;
		meta[1]=modelSize;
		meta[2]=st.foldCount;
		meta[3]=st.xvConst;

		strcpy(fn,filename);
		

		strcpy(label,st.label.c_str());
		
	
		strcpy(consistency,st.xvConsistency.c_str());
		log<<label<<" "<<consistency<<" "<<filename<<"\n";
		
		for (uint i=0; i<st.foldCount; i++) {
			log<<"\tFold: "<<i+1<<"\t"<<st.folds[i].label<<"\t";
			log<<setw(9)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<st.folds[i].clsError<<"%";
			log<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<st.folds[i].trainingFitness;
			log<<setw(9)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<st.folds[i].predError<<"%";
			log<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<st.folds[i].testingFitness;
			log<<setw(10)<<st.folds[i].oddsRatio<<"\n";
		}
		MPI::COMM_WORLD.Send(fn, 4096, MPI::CHAR, masterID, RESULTS_FILENAME);
		MPI::COMM_WORLD.Send(consistency, 24, MPI::CHAR, masterID, RESULTS_XV_CONSISTENCY);
		MPI::COMM_WORLD.Send(label, 128, MPI::CHAR, masterID, RESULTS_LABEL);
		
		MPI::COMM_WORLD.Send(stats, 5, MPI::FLOAT, masterID, RESULTS_STATS);
		MPI::COMM_WORLD.Send(stats, 5, MPI::INT, masterID, RESULTS_STATS);
		MPI::COMM_WORLD.Send(meta, 4, MPI::INT, masterID, RESULTS_META);
	}
#endif
	REPORTUNLOCK;
}






}

}
