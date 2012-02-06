//
// C++ Interface: richardslogistic
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_POPULATIONGROWTHRICHARDSLOGISTIC_H
#define SIMULATION_POPULATIONGROWTHRICHARDSLOGISTIC_H
#include <string>
#include "logisticgrowth.h"

namespace Simulation {
namespace PopulationGrowth {

using namespace std;

class RichardsLogistic : public LogisticGrowth {
public:
	~RichardsLogistic() {}

	//RichardsLogistic(uint initPopulationSize, uint carry);
	RichardsLogistic(uint initPopulationSize, uint carry, float m, float growthRate, float t, float var, uint reset = 0);
//	RichardsLogistic( uint initPopulationSize, uint carry, float , float growthRate) : 
		//LogisticGrowth(initPopulationSize, 0.0, carry-initPopulationSize, 0.0), timeOfMaxGrowth(0.0), polarity(0.0) { }

	uint GetBasePopulation(uint t);

	/**
	 * @brief Used for configuration summary
	 */
	void GenerateReport(std::ostream& os, uint headerWidth);

	/**
	 * @brief returns a string representation of the growth type
	 */
	string GetType() { return "Richard's Logistic"; }

	string GetHtmlReport();


	void SetTimeOfMaxGrowth(float tomg) { timeOfMaxGrowth = tomg; }
	void SetPolarity(float pol) { polarity = pol; }
	float GetTimeOfMaxGrowth() { return timeOfMaxGrowth; }
	float GetPolarity() { return polarity; }

	string GenerateCfgString();
	
protected:
	float timeOfMaxGrowth;						///<adjusts when we experience peak growth
	float polarity;								///<Affects which asymptote max growth occurs
	float resetPoint;							///<Which generation afterwhich a reset occurs (bottleneck)

};

/*** Richard's Logistic Growth **************************************/
inline
RichardsLogistic::RichardsLogistic( uint initPopulationSize, uint carry, float m, float growthRate, float t, float var, uint reset) : 
		LogisticGrowth(initPopulationSize, growthRate * -1, carry, var), timeOfMaxGrowth(m), polarity(t), resetPoint(reset) {
	modelType = RichardsLogisticGrowthRate;
}

inline
uint RichardsLogistic::GetBasePopulation(uint t) {
	if (t > resetPoint)
		t-=resetPoint;
	return (uint)(initPopulationSize + (carry / pow(1.0 + (polarity * exp(growthRate*(t-timeOfMaxGrowth))),(1.0/polarity))));
}

}

}

#endif
