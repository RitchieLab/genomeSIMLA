//
// C++ Implementation: snpsearchapplication
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snpsearchapplication.h"
#include <boost/timer.hpp>
#include "genetics/snprepostxtsorted.h"
#include "genetics/modelstatistics.h"
#include "evalbalancedaccuracypdt.h"


namespace MDR {

#ifdef USE_MPI
pthread_mutex_t SnpSearchApplication::input_lock 	= PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t SnpSearchApplication::output_lock 	= PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t SnpSearchApplication::mpi_lock		= PTHREAD_MUTEX_INITIALIZER;
#endif
using namespace boost;


SnpSearchApplication::SnpSearchApplication() : curTest(1), fileParser(NULL), fDist(NULL), oDist(NULL), pDist(NULL), nodeID(0), masterID(0), processCount(1)
	{
/*
#ifdef USE_MPI
	pthread_mutex_init(&input_lock, NULL);
	pthread_mutex_init(&output_lock, NULL);
	pthread_mutex_init(&mpi_lock, NULL);

#endif 	//USE_MPI
*/
	configuration = EseConfiguration::Instance();
}


SnpSearchApplication::~SnpSearchApplication()
{
}



/**	
 * Starts execution
 */
void SnpSearchApplication::MpiPrep()	{
#ifdef USE_MPI
	MPILOCK;
	nodeID = MPI::COMM_WORLD.Get_rank();
	processCount = MPI::COMM_WORLD.Get_size();
	MPIUNLOCK;
#endif

}

void SnpSearchApplication::MpiClose() {
#ifdef USE_MPI
	MPI::Finalize();
#endif
}



bool SnpSearchApplication::NextTestID(int &testID) {
	bool doContinue = false;
	TESTLOCK;
	if (nodeID == masterID) {
		if (curTest <= configuration->pTestCount)  {
			testID = curTest++;
			doContinue = true;
		}
	} else {
#ifdef USE_MPI
		MPI::COMM_WORLD.Recv(&testID, 1,  MPI::INT, masterID, TESTID);
		doContinue = testID!=0;
#else
		doContinue = false;
#endif

	}
	TESTUNLOCK;
	return doContinue;
}



float SnpSearchApplication::RunPTest(int testSeed, int &modelCount, vector<ModelStatistics> &stats) {
	//Let's make sure we are reporting only statistics for this test
	stats.clear();

	timer progress;
	SnpRepository ptest;

	//This is a quick fix for random numbers. 
	Utility::Random::globalGenerator.Seed(configuration->pTestSeed + testSeed);

	SnpReposTxtSorted bestScore(configuration->reportModelCount, configuration->crossValidationCount);

	//Acquire the evaluation suite- this time, we want randomized status
	SnpEvalSuite *pEval = configuration->GetPdtEval(true, testSeed);

	//If this fails, the application didn't setup the fileparser object (which must happen before this point)
	assert(fileParser);

	//Populate the test repository
	ptest.LoadData(fileParser);
	//Acquire the status. Technically, this should not change, but, in case the order changed for some reason, we want to make 
	//sure we are up to date.
	CaseControlStatus stat;
	fileParser->GetStatusMask(&stat);

	log<<setw(6)<<testSeed<<setw(15)<<setprecision(3)<<progress.elapsed();
	progress.restart();
	//Perform search on test data
	modelCount = ptest.Evaluate( configuration->comboStart - 1, configuration->comboEnd - 1, NULL, pEval);
	float timeToEvaluate = (float)progress.elapsed();
	
	bestScore.Sort();

	stats.clear();

	for (uint modelSize= 0; modelSize < configuration->comboEnd; modelSize++) {
		if (modelSize >= configuration->comboStart-1) {
			SnpAligned* bestModel = NULL;
			Genetics::Evaluation::ModelStatistics st = pEval->GetTopModel(modelSize, bestModel);
			st.GatherTestStatistics(&ptest, (EvalBalancedAccuracyPDT *)pEval->GetTrainer(), NULL);

			//cout<<"Node #"<<nodeID<<" ("<<testSeed<<") "<<modelSize<<" Model :"<<st.GetLabel()<<" "<<st.GetAvgTesting()<<"\n";

			stats.push_back(st);
		}
	}
	//pEval->UpdateDistribution();
	log<<setw(30)<<progress.elapsed()<<setw(15)<<modelCount<<"\n";
	delete pEval;
	
	return timeToEvaluate;
}

/**
 * I replaced the STL use of RAND to use the Utility::Random object, so this isn't necessary any more
 */
bool SnpSearchApplication::PerformSanityCheck() {

	return true;

}


}
