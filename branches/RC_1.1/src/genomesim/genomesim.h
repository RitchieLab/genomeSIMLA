//
// C++ Interface: genomesim
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIMGENOMESIM_H
#define GENOMESIMGENOMESIM_H

#include "utility/application.h"
#include "config.h"

namespace GenomeSIM {

/**
@brief Application to control the simulation

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenomeSim : public Utility::Application{
public:
    GenomeSim();

    ~GenomeSim();

	/**
	 * @brief Prints the help contents
	 */
	virtual void PrintHelp();
	
	/**	
	 * @brief Starts execution
	 */
	virtual void Start();

	/**
	 * @brief Parse the command line arguments
	 */
	bool ParseCmdLine(int argc, char **argv);

	/**
	 * @brief Parses the arguments past the configuration file
	 */
	int ParseCmd(int curr, int argc, char **argv);

	/**
	 * @brief Output the datasets to filesystem	
	 */
	void WriteDatasets();

protected:
	Config configuration;						///<The raw configuration
	uint startGeneration;						///<Used to determine if we start with a previous set
	uint additionalGenerations;					///<How many more generations are to be run
	bool GetBoolean(const char *value);			///<Converts value to a boolean value (On/off, Yes/No, etc)
	string OpenSummaryReport(ofstream &summary);
	
};

}

#endif
