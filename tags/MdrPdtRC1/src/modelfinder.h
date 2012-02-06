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
#include "appconfiguration.h"
#include "genetics/snprepository.h"
#include "genetics/snprepostxtsorted.h"
namespace MDRPDT {

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

	void PedigreeDescribeModels();

// 	
	bool DumpSampleConfig(int argc, char **argv);
	/** 
	 * @brief Parses the command line and sets up any local data
	 */
	bool ParseCmdLine(int argc, char **argv);
protected:
	bool VerifyConfiguration();					///<Do a quick sanity check on the settings provided
	string 					configFile;			///<Configuration file 
	AppConfiguration 	   	*configuration;		///<The raw configuration
	SnpRepository 			snps;				///<The original set of snps
	vector<string> 			models;				///<List of models to be described in detail
	bool 					doAnalyze;			///<Indicates whether to analyze or describe models
};

}

#endif
