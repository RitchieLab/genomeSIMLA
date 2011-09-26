//
// C++ Implementation: exponentialgrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "exponentialgrowth.h"

namespace Simulation {
namespace PopulationGrowth {

/*** Exponential Growth *************************************************/
ExponentialGrowth::ExponentialGrowth(uint initSize, float rate) : initPopSize(initSize), rate(rate) {
	r0=log(1+rate);
}

uint ExponentialGrowth::GetPopulationSize( uint generation ) {
	double newSize = AdjustedGrowth(pow(initPopSize, (r0 * (float)generation)) + (float)initPopSize);
	
	if (newSize > maxPoolSize)
		return maxPoolSize;
	else if (newSize < minPoolSize)
		return minPoolSize;
	else
		return (uint)newSize;
}

uint ExponentialGrowth::operator ( )( uint generation ) {
	return GetPopulationSize( generation );
}

uint ExponentialGrowth::GetInitialPopulationSize() {
	return initPopSize;
}

void ExponentialGrowth::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Population Growth Model: "<<"Exponential Growth\n";
	os<<setw(headerWidth)<<"Initial Population Size: "<<initPopSize<<endl;
	os<<setw(headerWidth)<<"Growth Rate: "<<setprecision(4)<<rate<<endl;
}


}

}
