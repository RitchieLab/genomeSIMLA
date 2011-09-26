//
// C++ Interface: growthmodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_GROWTHRATEGROWTHRATE_H
#define SIMULATION_GROWTHRATEGROWTHRATE_H

#include <math.h>
#include "utility/random.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace Simulation {

namespace PopulationGrowth {

enum GrowthRateType 	{	LinearGrowthRate = 0, ExponentialGrowthRate, LogisticGrowthRate, RichardsLogisticGrowthRate, GrowthRateUnknown };

/**
 * @brief Base class for the different possible growth rates
 * By default, we will make this model return the max size for each generation. It will 
 * effectively represent no growth at all and serve as the default setting
 @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GrowthRate{
public:
    GrowthRate(float variation);
    virtual ~GrowthRate();

	/**
	 * @brief Return the initial population size
	 */
	virtual uint GetInitialPopulationSize();

	/**
	 * @brief Returns the population size at generation (time) t
 	 */
	virtual uint GetPopulationSize(uint t);

	/**
	 * @brief Allow these to be used like a function
	 */
	virtual uint operator()(uint t);

	static uint maxPoolSize;				///<App wide max size. This is valid for ALL models
	static uint minPoolSize;				///<App wide min
//	static float variation;					///<Percentage of variation that is expected to occur during each generation
	static uint randomSeed;
	static uint adjustmentPeriod;
	/**
	 * Let's make the random number generation be easy to get/set. To start with, I'll just grab the 
	 * global one, but eventually it should be set up properly	
	 */
	Utility::Random rnd;

	/**
	 * @brief Quick adjustment to simulate natural variations in growth 
	 */
	double AdjustedGrowth(double population);

	/**
	 * @brief Used for configuration summary
	 */
	virtual void GenerateReport(std::ostream& os, uint headerWidth);

	virtual string GenerateCfgString() = 0;

	/**
	 * @brief Loads a model based on the contents of the input stream
	 */
	static GrowthRate *LoadModel(std::istream&os);

	/**
	 * @brief Dump an example of the growth the is expected from the local function. 
	 * @note These values do not represent the currently mimic the random draws that will be made. 
	 * As a result the resulting population sizes at a given generation should be close, they will most
	 * likely be at least somewhat different.
	 */
	void DiagramGrowth(std::ostream& os, int start, int stop, int interval);

	virtual string GetHtmlReport();

	/**
	 * @brief Returns a string representation of the growth type
	 */
	virtual string GetType();

	void SetVariation(float var);
	float GetVariation();

	string DrawGrowthChart(const char *project, int start, int stop, int interval, uint width = 800, uint height = 400, bool returnPage = true);

	GrowthRateType GetModelType();

protected:
	float variation;				///<Percentage of variation that is expected to occur during each generation
	float *draws;					///<Array of draws
	GrowthRateType modelType;
};

inline 
GrowthRateType GrowthRate::GetModelType() { 
	return modelType;
}








}

}

#endif
