//
// C++ Interface: logisticgrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_POPULATIONGROWTHLOGISTICGROWTH_H
#define SIMULATION_POPULATIONGROWTHLOGISTICGROWTH_H

#include "utility/types.h"
#include "growthrate.h"

namespace Simulation {

namespace PopulationGrowth {

/**
 * @brief Logistic Growth model (Principles of Population Genetics (hartl, clark) eq 1.12)
 * @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class LogisticGrowth : public GrowthRate {
public:
	/**
	 * @brief construction
	 * @param initPopulationSize The starting size of the population
	 * @param growthRate the growth rate 
	 * @param carry The Carrying capacity of the "environment"
	 */
	LogisticGrowth(uint initPopulationSize, float growthRate, uint carry, float var);
	~LogisticGrowth();

	/**
	 * @brief Returns the population size at generation (time) t
 	 */
	uint GetPopulationSize(uint t);
	virtual uint GetBasePopulation(uint t);
	/**
	 * @brief Allow these to be used like a function
	 */
	uint operator()(uint t);

	/**
	 * @brief Used for configuration summary
	 */
	void GenerateReport(std::ostream& os, uint headerWidth);

	/**
	 * @brief Returns the population size at time (0)
	 */
	uint GetInitialPopulationSize();

	/**
	 * @brief returns a string representation of the growth type
	 */
	string GetType() { return "Logistic"; }
	uint FindInflection(uint startingPoint, int maxTries);
	
	int GetCarryCapacity() { return (int)carry; }
	float GetGrowthRate() { return growthRate; }

	void SetCarryCapacity(float cap) { carry = cap; 	C=(carry - initPopulationSize)/initPopulationSize; }
	void SetInitialPopulationSize(float size) { initPopulationSize = size; 	C=(carry - size)/size;}
	void SetGrowthRate(float rate) { growthRate = rate; }

	string GenerateCfgString();

	string GetHtmlReport();
protected:
	float initPopulationSize;				///<The starting size of the population
	float growthRate;						///<Rate of growth
	float carry;							///<Carrying capacity of the environment
	float C;								///<This is a constant derived from the initial population
};




/*** Logistic Growth *************************************************/

inline
LogisticGrowth::LogisticGrowth(uint size, float rate, uint carry, float var) : 
		GrowthRate(var), initPopulationSize((float)size), growthRate(rate), carry((float)carry) { 
	modelType  = LogisticGrowthRate;

	if (initPopulationSize > maxPoolSize)
		initPopulationSize = maxPoolSize;
	C=(carry - initPopulationSize)/initPopulationSize;
}
inline
LogisticGrowth::~LogisticGrowth() {} 

inline
uint LogisticGrowth::GetBasePopulation(uint t) {
	return (uint)(carry/(1.0 + (C * exp(-1.0 * t * growthRate))));
}

inline
uint LogisticGrowth::GetPopulationSize( uint t) {
	//double basePop = carry/(1.0 + (C * exp (-1.0 * t * growthRate)));
	uint newSize = (uint)AdjustedGrowth(GetBasePopulation(t));
	if (newSize > maxPoolSize)
		newSize = maxPoolSize;
	else if (newSize < minPoolSize)
		newSize = minPoolSize;

	return newSize;
}

inline
uint LogisticGrowth::GetInitialPopulationSize() {
	return (uint)initPopulationSize;
}

inline
uint LogisticGrowth::operator ( )( uint t) {
	return GetPopulationSize(t);
}

}

}

#endif
