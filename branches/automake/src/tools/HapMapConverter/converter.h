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
#include "locusconverter.h"
#include "individual.h"
#include <map>

namespace Tools {

using namespace std;
/**
@brief Application to control the simulation

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Converter : public Utility::Application{
public:
    Converter();

    ~Converter();

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
	 * @brief Output the datasets to filesystem	
	 */
	void WriteData();

	void PurgeLoci();
	void PurgeIndividuals();

protected:
	Config configuration;						///<The raw configuration

	vector<LocusConverter *> loci;				///<Each of the loci we will be converting
	vector<Individual *>individuals;			///<All of the people of interest
	map<string, LocusConverter*> locusLookup;	///<In order to flip the doWrite bits on and off


	/**
	 * @brief Load the legend data from text file
	 */
	void LoadLegend();
	void LoadHaplotypes();
	void LoadSampleDetails();
	void ParseInclusionsList();
	
};

}

#endif
