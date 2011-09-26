//
// C++ Interface: eseapp
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEESEAPP_H
#define ESEESEAPP_H
#include "utility/utility.h"
#include "eseconfiguration.h"
#include "genetics/snprepository.h"
#include "genetics/snprepostxtsorted.h"


#define TESTID 1		//This represents the message ID for test ids
#define RESULT 2		//This represents the message ID for results
//		RESULT     		(INT)    	This is the test ID
//		RESULT[0]  		(FLOAT)      -- Loadtime for dataset
//      RESULT[1]  		(FLOAT)      -- Runtime for test
//      RESULT[3..N+2] 	(FLOAT)  	This is the array of scores
//      RESULT 			(CHAR) 		This is the set of model IDs

namespace MDR {

using namespace Utility;
using namespace Genetics;

/**
@brief The main application for ModelFinder

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ModelFinder : public Application {
public:
    ModelFinder();
    ~ModelFinder();

	/**
	 * @brief Prints the help contents
	 */
	virtual void PrintHelp();
	
	/**	
	 * @brief Starts execution
	 */
	virtual void Start();

	/**
	 * @brief Perform the pedigree search
	 */
	virtual void PedigreeSearch();
	
	void PedigreeSlaveSearch();

	void PedigreeDescribeModels();

	/**
	 * @brief Perform a case/control based search
	 */
	virtual void CaseControlSearch();
	
	bool DumpSampleConfig(int argc, char **argv);
	/** 
	 * @brief Parses the command line and sets up any local data
	 */
	bool ParseCmdLine(int argc, char **argv);

	bool NextTestID(int &testID);

	//void AppendPTest(int testID, int modelCount, float *pResult, PTestDistribution *dist, const char *modelIDs);



#ifdef USE_MPI
	static void *ProcessWatcher(void *arg);
#endif
protected:
	bool VerifyConfiguration();					///<Do a quick sanity check on the settings provided
	string 					configFile;			///<Configuration file 
	EseConfiguration 	   	*configuration;		///<The raw configuration
	SnpRepository 			snps;				///<The original set of snps
	vector<string> 			models;				///<List of models to be described in detail
	bool 					doAnalyze;			///<Indicates whether to analyze or describe models

	uint 					curTest;			///<The test currently in use


	//MPI stuff
	int 					nodeID;				///<Used for MPI self ID
	int 					masterID;			///<Used for Master ID
	int						processCount;		///<Count of processes associated with this job
	
	
	
};

}

#endif
