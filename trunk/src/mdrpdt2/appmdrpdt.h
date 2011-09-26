//
// C++ Interface: appmdrpdt
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTAPPMDRPDT_H
#define MDRPDTAPPMDRPDT_H
#include "tstatistic.h"
#include "ptestdistribution.h"
#include "familyrepository.h"
#include "genotyperepository.h"
#include "appconfig.h"
#include <pthread.h>

#ifdef USE_MPI
#include "mpi.h"

#define TESTLOCK pthread_mutex_lock(&input_lock)
#define TESTUNLOCK pthread_mutex_unlock(&input_lock)
#define REPORTLOCK pthread_mutex_lock(&output_lock)
#define REPORTUNLOCK pthread_mutex_unlock(&output_lock)
#define MPILOCK pthread_mutex_lock(&mpi_lock)
#define MPIUNLOCK pthread_mutex_unlock(&mpi_lock)

#else

#define TESTLOCK
#define TESTUNLOCK
#define REPORTLOCK
#define REPORTUNLOCK
#define MPILOCK
#define MPIUNLOCK

#endif

/**
 * MPI MESSAGES
 * 
 * MPI_TESTID
 * 		TestID    	INT(1)			0 indicates stop processing
 * 
 * MPI_RESULT
 * 		TestID		INT(1)
 * 		ModelID		CHAR(MAX_SIZE)	
 * 		RESULT_MOR	FLOAT(1)			For now, let's just send the value of the best. If we need all N scores, we'll have to change things a bit
 */

namespace MdrPDT {


#define MPI_TESTID 		1
#define MPI_SHUTDOWN 	2
#define MPI_RESULT 		3
/**
@brief The main application object

	@author Eric Torstenson
*/
class AppMdrPDT{
public:


    AppMdrPDT();
    ~AppMdrPDT();

	/**
	 * @brief Load the pedigree data. This must be done prior to any evaluation
	 */
	void LoadData(const char *dataset);

	/**
	 * @brief Evaluates a given dataset using the method, eval, and stores results into report
	 */
	void EvaluateDataset(Evaluation::TFinalReport& report);

	/**
	 * @brief debug the PTest
	 */
	void EvaluatePTest(int testNumber, Evaluation::TFinalReport& report);

	/**
	 * @brief Run the ptests, returning the resulting distribution
	 */
	Distribution::PTestDistribution *RunPTests(int ptestCount);

	/**
	 * @brief Return the evaluation method, in case it is needed for reporting
	 */
	Evaluation::EvaluationMethod *GetEvaluationMethod();

	/**
	 * @brief Loads the configuration from the file, cfgFilenme
	 */
	ApplicationConfiguration *LoadConfiguration(const char *cfgFilename);

	/**
	 * @brief Return the seed in use
	 */
	int GetSeed() { return seed; }

	/**
	 * @brief Set the seed for permutation testing
	 */
	void SetSeed(int seed);

	/**
	 * @brief Return the number of folds to be used in analysis (Cross validation)
	 */
	int GetFoldCount() { return xvCount; }
	/**
	 * @brief Pass the arguments to the application object
	 */
	bool ParseCmdLine(int argc, char **argv);
	void PrintHelp();				///<Display usage details
	void PrintBanner();				///<Display details about the software

	void SetBestModel(float bestMOR, const char *bestModel, int ptestCount);
	/**
	 * @Abstract out the method for filename creation. We might switch to setting these up differently later
	 */
	std::string GetFilename(const char *extension);

	/**
	 * @brief returns next seed, or -1 if no more tests are to be done
	 */
	int GetNextPTest();
	void ResetPTestCounts(int ptestCount);

	bool IsMasterNode() { return masterID == nodeID; }
#ifdef USE_MPI
	void AppendToMORDist(Distribution::PTestDistribution* disto, int testID, const char *model, float score);
	static void *ProcessWatcher(void*);
#endif
protected:
	int seed;						///<Random number seed
	int xvCount;					///<The number of cross folds required for XV
	string targetModel;				///<The model to be described (if not typical analysis)
	/**
	 * @brief The pedigree repository
	 */
	FamilyRepository<MdrPDT::GenotypeConversion> repo;

	/**
	 * @brief This is the raw data used for evaluating the original data
	 */
	GenotypeRepository data;				

	Evaluation::EvaluationMethod *eval;			///<Evaluation object (T-statistic)
	
	ConfigurationReader cfg;					///<Configuration parser/recording

	ofstream pedigreeReport;					///<Pedigree report stream
	ofstream mpiReport;							///<Each node will have it's own log

	float winningScore;							///<Score for short circuiting
	string winningID;							///<Best model ID (for short circuiting report)
	int ptestShortCircuit;						///<Number of failed tests required for short circuiting

	int 					ptestCount;			///<The number of permutation tests
	int 					ptestsRequested;	///<Keep count of how many we've done
	int 					ptestExceptions;	///<Used to count of how many tests have exceeded the winning model
	int 					threadCount;		///<Number of independant threads to use

	static void *RunPTests(void *args);
	void RunPTests(Distribution::PTestDistribution*);

	int 					nodeID;				///<Used for MPI self ID
	int 					masterID;			///<Used for Master ID
	int						processCount;		///<Count of processes associated with this job

	void MpiPrep();
	void MpiClose();
#ifdef USE_MPI
	static pthread_mutex_t input_lock;			///<Used to lock ptest input
	static pthread_mutex_t output_lock;			///<Used to lock ptest ouput
	static pthread_mutex_t mpi_lock;			///<Used to lock MPI calls

	std::streambuf* origCout;					///<Used to cache the cout for slave nodes
#endif

};


}

#endif
