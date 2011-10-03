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
	os<<setw(headerWidth)<<"Variation: "<<setprecision(4)<<variation<<endl;
}

string LinearGrowth::GetHtmlReport() {
	stringstream ss;
	ss<<"<TABLE><TR><TH>Initial Population</TH><TH>Growth Rate</TH><TH>Variation></TH><TH>Min. Population</TH><TH>Max. Population</TH></TR>\n";
	ss<<"\t<TR><TD>"<<initPopulationSize<<"</TD><TD>"<<growthRate<<"</TD><TD>"<<variation<<"</TD><TD>"<<minPoolSize<<"</TD><TD>"<<maxPoolSize<<"</TD></TR></TABLE>\n";
	return ss.str();
}

string LinearGrowth::GenerateCfgString() {
	stringstream ss;
	ss<<"GROWTH_RATE LINEAR "<<initPopulationSize<<" "<<variation<<" "<<growthRate;
	return ss.str();
}

}

}
