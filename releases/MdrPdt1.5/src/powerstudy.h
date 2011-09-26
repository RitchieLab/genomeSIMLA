//
// C++ Interface: powerstudy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPOWERSTUDY_H
#define MDRPOWERSTUDY_H
#include "utility/utility.h"
#include "eseconfiguration.h"
#include "genetics/snprepository.h"


namespace MDR {

namespace Power {
/**
@brief Application designed to perform power studies on MDR-PDT

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PowerStudy : public Application {
public:

struct ModelResults {
	ModelResults() : freq(0), totalT(0.0) { }

	double GetAverageT() {
		return totalT / (float)freq;
	}

	uint freq;
	float totalT;
	
};

    PowerStudy();

    ~PowerStudy();



	/**
 	 * @brief This will let the program set up the necessary parameters based on cmd line arguments. 
	 * @return true indicates that execution can begin
 	 */
	virtual bool ParseCmdLine(int argc, char** argv);

	/**
	 * @brief Prints the help contents
	 */
	virtual void PrintHelp();
	
	/**	
	 * @brief Starts execution
	 */
	virtual void Start();
	
	string NextDataset();

protected:

	/**
	 * @brief Do a quick sanity check on the settings provided
	 */
	bool VerifyConfiguration();
	string targetModel;
	vector<string> inputfiles;

	map<string, ModelResults> power;
	uint fileIdx;								///<Used to keep up with the index of the current file
	string configFile;
	EseConfiguration 	   	*configuration;		///<The raw configuration
	double PedigreeSearch(const char *filename);

};

}
}
#endif
