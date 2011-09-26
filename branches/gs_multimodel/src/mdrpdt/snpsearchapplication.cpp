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

string SnpSearchApplication::sasDirectory = "sasdir";
string SnpSearchApplication::sasApp		  = "~/bin/sas";
bool SnpSearchApplication::generateRegressionScript = false;

using namespace boost;


SnpSearchApplication::SnpSearchApplication() : curTest(1), fileParser(NULL), fDist(NULL), oDist(NULL), pDist(NULL), lrDist(NULL), nodeID(0), masterID(0), processCount(1), DoWriteSAS(false)
	{
/*
#ifdef USE_MPI
	pthread_mutex_init(&input_lock, NULL);
	pthread_mutex_init(&output_lock, NULL);
	pthread_mutex_init(&mpi_lock, NULL);

#endif 	//USE_MPI


*/
	configuration = EseConfiguration::Instance();
	pTestCount = 0;
	seeds = NULL;
	rand = NULL;
}


void SnpSearchApplication::InitRandom() {
	pTestCount = configuration->pTestCount;
	seeds = new int[pTestCount + 1];

	rand = &Utility::Random::globalGenerator;
	rand->Seed(configuration->pTestSeed);

	assert(pTestCount >= 0);
	//Set up the random numbers
	for (size_t i=0; i<(size_t)pTestCount; i++) {
		seeds[i]=rand->lrand();
//		cout<<"seeds["<<i<<"] = "<<seeds[i]<<"\n";
	}
}

SnpSearchApplication::~SnpSearchApplication()
{
	delete[] seeds;
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



bool SnpSearchApplication::NextTestID(size_t &testID) {
	bool doContinue = false;
	TESTLOCK;
		
	if (seeds == NULL) 
		InitRandom();
	if (nodeID == masterID) {
		if (curTest <= (size_t)pTestCount)  {
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




float SnpSearchApplication::CalculateLikelihood(const char *filename) {
	float likelihood = 0.0;

	string reportFilename = string(filename) + ".report";
	char cmd[2048];
	sprintf(cmd, "%s %s -log %s.log -print %s", sasApp.c_str(), filename, filename, reportFilename.c_str());
	int rv = system(cmd);
	
	if (rv != 0) {
		cout<<"There was an error when trying to run sas\n";
		return -1.0;
	}
	
	ifstream file(reportFilename.c_str(), ios_base::in);
	char line[2048];
	string junk;
	
	for (size_t i=0; i<27; i++)
		file.getline(line, 2048);

	file>>junk>>junk>>junk>>junk>>likelihood;

	cout<<"Found the following in the report: "<<likelihood<<"\n";

	return likelihood;
}

float SnpSearchApplication::WriteSasScripts(SnpRepository &snps, SnpEvalSuite *e, SnpAligned *bestModel, const char *filename, const char *type) {
	stringstream saScript;						///<Name of the sas script
	stringstream sasDataFile;					///<Name of the specilized ped file
	stringstream sasIntScript;					///<Name of the interaction script
	string topmodel;

	
	//Make it known we are using it
	bestModel->IncrementInstanceCount();

	uint modelSize = bestModel->GetLabelCount();

	sasDataFile<<sasDirectory<<"/"<<filename<<"."<<nodeID<<topmodel<<"."<<modelSize<<"L.sasped";
	saScript<<sasDirectory<<"/"<<filename<<topmodel<<"."<<modelSize<<"L.sas";
	sasIntScript<<sasDirectory<<"/"<<filename<<topmodel<<"."<<modelSize<<"L.int.sas";

	regressionlog<<filename<<"("<<bestModel->GetLabel()<<")\t"<<type<<"\t"<<saScript.str()<<"\t"<<sasIntScript.str()<<"\n";


	vector<SnpAligned *> composites;

	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)e->GetTrainer();
	//The line parser is needed for writing the selected model snps back to file
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();


	for (uint snpID=0; snpID < modelSize; snpID++) {
		char lbl[128];
		uint pos = bestModel->GetLabel(snpID);
		sprintf(lbl, "%d", pos);
		SnpAligned *m = snps.GetSnp(lbl);
		//We want to make sure these have data from these if we didn't do single runs

		evalPdt->GetHighRiskCells(m);

		composites.push_back(snps.GetSnp(lbl));
	}

if (DoWriteSAS) {

//	cout<<"Getting ready to write sas datafile: "<<sasDataFile.str()<<"\n";

	//Write the minimal ped file
	pdtParser->WriteModelLociOnly(bestModel->GetLabel(), sasDataFile.str().c_str());
	//I Think I will need to reconstitute the model for this part
	
	pdtParser->WriteSaScript(saScript.str().c_str(), composites, sasDataFile.str().c_str(), false);

	pdtParser->WriteSaScript(sasIntScript.str().c_str(), composites, sasDataFile.str().c_str(), true);
}
	float likelihoodRatio = 0.0;			//= CalculateLikelihood(sasIntScript.str().c_str()) - CalculateLikelihood(sasDataFile.str().c_str());

	for (uint snpID=0; snpID<composites.size(); snpID++) 
		composites[snpID]->ReduceInstanceCount();

	//Make sure we don't leave this dangling
	bestModel->ReduceInstanceCount();
	
	return likelihoodRatio;
}


float SnpSearchApplication::RunPTest(int testSeed, int &modelCount, vector<ModelStatistics> &stats) {
	//Let's make sure we are reporting only statistics for this test
	stats.clear();

	timer progress;
	SnpRepository ptest;

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

//			cout<<"Node #"<<nodeID<<" ("<<testSeed<<") "<<modelSize<<" Model :"<<st.GetLabel()<<" "<<st.GetAvgTesting()<<"\n";

			stats.push_back(st);
			if (bestModel)
				bestModel->ReduceInstanceCount();
		}
		else
			stats.push_back(Genetics::Evaluation::ModelStatistics());
	}



/**
 * Here we are setting up the SAS ptests
 */
	ModelStatistics &bestMOR = stats[SelectBestMOR(stats)];
	ModelStatistics &bestPE  = stats[SelectBestPE(stats)];

	char filename[2048];

/*
//Quickly verify that the ptests are different
	sprintf(filename, "ptest-%d.ped", testSeed);
	EvalBalancedAccuracyPDT *evalPdt = (EvalBalancedAccuracyPDT *)pEval->GetTrainer();
	GtLineParserMdrPdt *pdtParser = (GtLineParserMdrPdt*)((GtFileParserBuffered *)fileParser)->GetLineParser();
	pdtParser->WritePedFile(filename);
//end debug verification
*/
	if (generateRegressionScript) {
		sprintf(filename, "ptest-%d-MOR", testSeed);
		SnpAligned *bestModel = ptest.GetSnp(bestMOR.label.c_str());
		if (bestModel) {
			bestMOR.likelihoodRatio = WriteSasScripts(ptest, pEval, bestModel, filename, "MOR-PTest");
			bestModel->ReduceInstanceCount();
		}
	
		sprintf(filename, "ptest-%d-PE", testSeed);
		bestModel = ptest.GetSnp(bestPE.label.c_str());
		if (bestModel) {
			bestPE.likelihoodRatio = WriteSasScripts(ptest, pEval, bestModel, filename, "PE-PTest");
			bestModel->ReduceInstanceCount();
		}
	}
	
/**
 * End SAS
 */

	//pEval->UpdateDistribution();
	log<<setw(30)<<progress.elapsed()<<setw(15)<<modelCount<<"\n";
	delete pEval;
	
	return timeToEvaluate;
}


/**
 * Scans through the different orders and determines which model is best based on
 * 		- T. Statistic & Cross validation frequency
 *      - Next (if tie) Matched Odds Ratio
 * 		- Returns the index to the winner
 */
size_t SnpSearchApplication::SelectBestMOR(vector<ModelStatistics>& stats) {
	string topORModel;

	stringstream tempReport;
	
	float bestOR = 100.0;
	int bestXV = 0;

	size_t topOR = (size_t)-1;
	uint orModelSize = 0;
	
//	uint statsCount = stats.size();

//	tempReport<<" --";
	//This just grabs the top model. 
	for (uint modelSize= configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
		if (modelSize < stats.size() ) {
			ModelStatistics &st = stats[modelSize];
		
			float curMOR = st.GetOddsRatio();
	
			if (modelSize == configuration->comboStart-1 || bestXV<=st.xvConst) {

//				tempReport<<" "<<st.GetLabel()<<" ["<<modelSize + 1<<"] MOR="<<curMOR<<" XVF="<<st.xvConst;
				//If we have better consistency than previous ones, set everything
				if (bestXV < st.xvConst) {
					bestOR = curMOR;
					topOR = modelSize;
					bestXV = st.xvConst;
					orModelSize = modelSize + 1;
				}
				else if (bestXV == st.xvConst) {
					if (curMOR > bestOR) {
						bestOR = curMOR;
						topOR = modelSize;
						orModelSize = modelSize + 1;
					}
	
				}
			}
		} else {
			cout<<"Hm, there are not enough models here! "<<modelSize<<" !< "<<stats.size()<<" \n";
			assert(modelSize < stats.size());
		}
	}
//	tempReport<<" = "<<topOR<<"/"<<bestXV;
//	cout<<tempReport.str()<<"\n";

	return topOR;
}

/**
 * Scans through the different orders and determines which model is best based on
 * 		- T. Statistic & Cross validation frequency
 *      - Next (if tie) Matched Odds Ratio
 * 		- Returns the index to the winner
 */
size_t SnpSearchApplication::SelectBestPE(vector<ModelStatistics>& stats) {
	string topORModel;
	string topPEModel;
	
	float bestPE = 100.0;
	int bestXV = 0;

	size_t topPE = (size_t)-1;
//	uint orModelSize = 0;
	
	stringstream tempReport;
//	tempReport<<" ++";

	//This just grabs the top model. 
	for (uint modelSize= configuration->comboStart-1; modelSize < configuration->comboEnd; modelSize++) {
		ModelStatistics &st = stats[modelSize];


		float curPE = st.GetAvgPredictionError();

//		tempReport<<" "<<st.GetLabel()<<" ["<<modelSize + 1<<"] PE="<<curPE<<" XVF="<<st.xvConst;

		if (modelSize == configuration->comboStart-1 || bestXV<=st.xvConst) {
			//If we have better consistency than previous ones, set everything
			if (bestXV < st.xvConst) {
				bestPE = curPE;
				topPE = modelSize;
				bestXV = st.xvConst;
			}
			else if (bestXV == st.xvConst && curPE < bestPE) {
				bestPE = curPE;
				topPE = modelSize;
			}
		}
	}
//	tempReport<<" = "<<topPE<<"/"<<bestXV;
//	cout<<tempReport.str()<<"\n";

	return topPE;
}

/**
 * I replaced the STL use of RAND to use the Utility::Random object, so this isn't necessary any more
 */
bool SnpSearchApplication::PerformSanityCheck() {

	return true;

}


}
