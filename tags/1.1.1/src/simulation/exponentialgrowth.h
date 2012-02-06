//
// C++ Interface: exponentialgrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPOPULATIONGROWTHEXPONENTIALGROWTH_H
#define SIMULATIONPOPULATIONGROWTHEXPONENTIALGROWTH_H
#include "growthrate.h"


namespace Simulation {
namespace PopulationGrowth {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
/**
 * @brief Exponential Growth model
 * @note Based on the formula from Principles of Population Genetics (hartl, clark) eq 1.6
 */
class ExponentialGrowth : public GrowthRate {
public:
	ExponentialGrowth(uint initPopulationSize, float growthRate, float var);
	~ExponentialGrowth() {}
	/**
	 * @brief Returns the population size at generation (time) t
 	 */
	uint GetPopulationSize(uint t);

	/**
	 * @brief Allow these to be used like a function
	 */
	uint operator()(uint t);

	/**
	 * @brief Returns the size of the population at T0
	 */
	uint GetInitialPopulationSize();

	/**
	 * @brief Used for configuration summary
	 */
	void GenerateReport(std::ostream& os, uint headerWidth);

	/**
	 * @brief Returns the textual description of the growth type
	 */
	string GetType() { return "Exponential"; }

	string GenerateCfgString();
	
	virtual string GetHtmlReport();

	float GetGrowthRate() { return rate; }
	void SetGrowthRate(float rate) { this->rate = rate; r0=log(1+rate);}
	void SetInitialPopulationSize(uint size) { initPopSize = size; }

protected:
	uint initPopSize;						///<N0
	float rate;								///<Variable growth rates based on min/max and random draw
	float maxGrowthRate;					///<Variable growth rates based on min/max and random draw
	double r0;								///<This is a pseudo constant based on the r and initial population size
};




}

}
#endif
