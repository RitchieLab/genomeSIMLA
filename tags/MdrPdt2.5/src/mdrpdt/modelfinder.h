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
#include "snpsearchapplication.h"
#include "genetics/snprepository.h"
#include "genetics/snprepostxtsorted.h"


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
class ModelFinder : public SnpSearchApplication {
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


	//void AppendPTest(int testID, int modelCount, float *pResult, PTestDistribution *dist, const char *modelIDs);
//	bool PerformSanityCheck();

#ifdef USE_MPI
	static void *ProcessWatcher(void *arg);
#endif
protected:
	bool VerifyConfiguration();					///<Do a quick sanity check on the settings provided
	string 					configFile;			///<Configuration file 
	SnpRepository 			snps;				///<The original set of snps
	vector<string> 			models;				///<List of models to be described in detail
	bool 					doAnalyze;			///<Indicates whether to analyze or describe models



	
	
	
};


}

#endif
