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

	bool NextTestID(int &testID);
	float RunPTest(int testSeed, int &modelCount, vector<ModelStatistics> &stats);
	bool PerformSanityCheck();
protected:
	EseConfiguration 	   	*configuration;		///<The raw configuration
	uint 					curTest;			///<The test currently in use
	GtFileParser 			*fileParser;		///<We need a file parser to be setup properly

	/**
	 * Copies of the distributions. 
	 */
	PTestDistribution 		*fDist;
	PTestDistribution 		*oDist;
	PTestDistribution 		*pDist;


	//MPI stuff
	int 					nodeID;				///<Used for MPI self ID
	int 					masterID;			///<Used for Master ID
	int						processCount;		///<Count of processes associated with this job
	std::ofstream 			log;


	void MpiPrep();
	void MpiClose();

	static pthread_mutex_t input_lock;
	static pthread_mutex_t output_lock;
	static pthread_mutex_t mpi_lock;


};

}

#endif
