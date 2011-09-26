//
// C++ Interface: snpsearchapplication
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRSNPSEARCHAPPLICATION_H
#define MDRSNPSEARCHAPPLICATION_H

#ifdef USE_MPI

#include "mpi.h"
#include <pthread.h>


#define TESTLOCK pthread_mutex_lock(&input_lock)
#define TESTUNLOCK pthread_mutex_unlock(&input_lock)
#define REPORTLOCK pthread_mutex_lock(&output_lock)
#define REPORTUNLOCK pthread_mutex_unlock(&output_lock)
#define MPILOCK pthread_mutex_lock(&mpi_lock)
#define MPIUNLOCK pthread_mutex_unlock(&mpi_lock)

#else
//We don't want these to mean anything if we arne't using MPI
#define TESTLOCK
#define TESTUNLOCK
#define REPORTLOCK
#define REPORTUNLOCK
#define MPILOCK
#define MPIUNLOCK
#endif 	//USE_MPI
#include <fstream>
#include <iostream>
#include "utility/application.h"
#include "genetics/ptestdistribution.h"
#include "eseconfiguration.h"


#define TESTID 1		//This represents the message ID for test ids
#define RESULT 2		//This represents the message ID for results
#define DATASET 3 		//This represents the message ID for a dataset


namespace MDR {


/**
This class should house the common functionality between the different search methods

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpSearchApplication : public Utility::Application {
public:
    SnpSearchApplication();

    ~SnpSearchApplication();

	/**
	 * @brief Sets the next test ID if there is one
	 * @Param testID Where the id will be stored
	 * @return Indicate that a test was returned (otherwise, we are finished)
	 */
	bool NextTestID(size_t &testID);
	
	/**
	 * @brief Runs a ptest 
	 * @param testSeed the seed to be used with the run
	 * @param modelCount Storage for the size of winning model in stats
	 * @param stats The statistics associated with the winner
	 * @return 
	 */
	float RunPTest(int testSeed, int &modelCount, vector<ModelStatistics> &stats);

	void InitRandom();
	/**
	 * @brief checks for any unreasonable conditions that should prevent execution
	 */
	bool PerformSanityCheck();

	/**
	 * @brief Writes (and executes) sas scripts and returns the likelihood ratio
	 * @param snps the repository where snps can be found
	 * @param e the evaluation object that is to be used
	 * @param bestModel The model to be used for calculating the likelihood ratio
	 * @param filename The name of the pedigree file
	 * @param isTopmodel (currently, this is a bit outdated)
	 * @return likelihood ratio (or 0.0)
	 */
	float WriteSasScripts(SnpRepository &snps, SnpEvalSuite *e, SnpAligned *bestModel, const char *filename, const char *type);
	float CalculateLikelihood(const char *filename);

	static bool generateRegressionScript;
	static string sasDirectory;
	static string sasApp;
protected:
	EseConfiguration 	   	*configuration;		///<The raw configuration
	size_t 					curTest;			///<The test currently in use
	GtFileParser 			*fileParser;		///<We need a file parser to be setup properly
	/**
	 * Copies of the distributions. 
	 */
	PTestDistribution 		*fDist;
	PTestDistribution 		*oDist;
	PTestDistribution 		*pDist;
	PTestDistribution	 	*lrDist;


	//MPI stuff
	int 					nodeID;				///<Used for MPI self ID
	int 					masterID;			///<Used for Master ID
	int						processCount;		///<Count of processes associated with this job
	std::ofstream 			log;
	ofstream regressionlog;
	bool DoWriteSAS;

	int *seeds;
	Utility::Random *rand;
	int pTestCount;

	void MpiPrep();
	void MpiClose();


#ifdef USE_MPI
	static pthread_mutex_t input_lock;			///<Used to lock ptest input
	static pthread_mutex_t output_lock;			///<Used to lock ptest ouput
	static pthread_mutex_t mpi_lock;			///<Used to lock MPI calls
#endif
	/**
 	 * @brief Selects the best model from stats based on PE as the tie breaker
	 */
	size_t SelectBestPE(vector<ModelStatistics>& stats);
	
	/**
	 * @brief Selects the best model from stats based on MOR as the tie breaker
	 */
	size_t SelectBestMOR(vector<ModelStatistics>& stats);
};

}

#endif
