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
	LinearGrowth(uint initSize, float growthRate);
	~LinearGrowth();
	
	/**
	 * @brief Returns the population size at generation (time) t
 	 */
	uint GetPopulationSize(uint t);

	/**
	 * @brief Return the initial population
	 */
	uint GetInitialPopulationSize();

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

protected:
	float initPopulationSize;				///<The starting size of the population
	float growthRate;						///<Rate of growth
};
/*** Linear Growth *************************************************/

inline
LinearGrowth::LinearGrowth(uint initSize, float rate) : initPopulationSize(initSize), growthRate(rate) { }

inline
LinearGrowth::~LinearGrowth()  { }

inline
uint LinearGrowth::GetPopulationSize( uint t) {
	return (uint)AdjustedGrowth(initPopulationSize + (growthRate * initPopulationSize * t));
}

uint LinearGrowth::GetInitialPopulationSize() {
	return (uint)initPopulationSize;
}


inline
uint LinearGrowth::operator ( )( uint t) {
	return GetPopulationSize(t);
}


}

}

#endif
