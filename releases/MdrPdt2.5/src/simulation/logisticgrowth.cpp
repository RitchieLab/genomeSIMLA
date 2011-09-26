//
// C++ Implementation: logisticgrowth
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "logisticgrowth.h"

namespace Simulation {

namespace PopulationGrowth {



void LogisticGrowth::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Population Growth Model: "<<"Logistic Growth\n";
	os<<setw(headerWidth)<<"Initial Population Size: "<<initPopulationSize<<endl;
	os<<setw(headerWidth)<<"Growth Rate: "<<setprecision(4)<<growthRate<<endl;
	os<<setw(headerWidth)<<"Carrying Capacity: "<<carry<<endl;
}




}

}
