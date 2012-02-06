//
// C++ Implementation: lineargrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineargrowth.h"

namespace Simulation {

namespace PopulationGrowth {


void LinearGrowth::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Population Growth Model: "<<"Linear Growth\n";
	os<<setw(headerWidth)<<"Initial Population Size: "<<initPopulationSize<<endl;
	os<<setw(headerWidth)<<"Growth Rate: "<<setprecision(4)<<growthRate<<endl;
}

}

}
