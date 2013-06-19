//
// C++ Interface: lineargrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_POPULATIONGROWTHLINEARGROWTH_H
#define SIMULATION_POPULATIONGROWTHLINEARGROWTH_H
#include "growthrate.h"

namespace Simulation {

namespace PopulationGrowth {

/**
 * @brief Simple linear growth
 * @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class LinearGrowth : public GrowthRate {
public:
	LinearGrowth(uint initSize, float growthRate, float var);
	~LinearGrowth();
	
	/**
	 * @brief Returns the population size at generation (time) t
 	 */
	uint GetPopulationSize(uint t);

	/**
	 * @brief Return the initial population
	 */
	uint GetInitialPopulationSize();

	void SetInitialPopulationSize(uint size);
	
	void SetGrowthRate(float gr);

	float GetGrowthRate() { return growthRate; }
	/**
	 * @brief Allow these to be used like a function
	 */
	uint operator()(uint t);

	/**
	 * @brief Used for configuration summary
	 */
	void GenerateReport(std::ostream& os, uint headerWidth);

	/**
	 * @brief Returns a string representation of the growth type
	 */
	string GetType() { return "Linear"; }

	/**
	 * @brief Generate the html report
	 */
	string GetHtmlReport();

	string GenerateCfgString();

protected:
	float initPopulationSize;				///<The starting size of the population
	float growthRate;						///<Rate of growth
};
/*** Linear Growth *************************************************/

inline
LinearGrowth::LinearGrowth(uint initSize, float rate, float var) : GrowthRate(var), initPopulationSize(initSize), growthRate(rate) { 	modelType = LinearGrowthRate;

	if (initPopulationSize > maxPoolSize)
		initPopulationSize = maxPoolSize;
}



inline
LinearGrowth::~LinearGrowth()  { }

inline
uint LinearGrowth::GetPopulationSize( uint t) {
	return (uint)AdjustedGrowth(initPopulationSize + (growthRate * t));
}

inline
uint LinearGrowth::GetInitialPopulationSize() {
	return (uint)initPopulationSize;
}

inline
void LinearGrowth::SetInitialPopulationSize(uint popSize) {
	initPopulationSize = popSize;
}

inline
void LinearGrowth::SetGrowthRate(float growthRate) {
	this->growthRate = growthRate;
}


inline
uint LinearGrowth::operator ( )( uint t) {
	uint newSize = GetPopulationSize(t);
	if (newSize > maxPoolSize)
		newSize = maxPoolSize;
	else if (newSize < minPoolSize)
		newSize = minPoolSize;
	return newSize;
}


}

}

#endif
