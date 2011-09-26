//
// C++ Implementation: richardslogistic
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "richardslogistic.h"
#include <sstream>

namespace Simulation {

namespace PopulationGrowth {

void RichardsLogistic::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Population Growth Model: "<<"Richard's Logistic Growth\n";
	os<<setw(headerWidth)<<"Initial Population Size: "<<initPopulationSize<<endl;
	os<<setw(headerWidth)<<"Growth Rate: "<<setprecision(4)<<((-1.0) * growthRate)<<endl;
	os<<setw(headerWidth)<<"Carrying Capacity: "<<(uint)(carry + initPopulationSize)<<endl;
	os<<setw(headerWidth)<<"Time of Max Growth: "<<timeOfMaxGrowth<<endl;
	os<<setw(headerWidth)<<"Polarity: "<<polarity<<endl;
	os<<setw(headerWidth)<<"Variation: "<<setprecision(4)<<variation<<endl;
}

string RichardsLogistic::GetHtmlReport() {
	stringstream ss;
	ss<<"<TABLE><TR><TH>Initial Population</TH><TH>Growth Rate</TH><TH>Carrying Capacity</TH><TH>Time of Max. Growth</TH><TH>Polarity</TH><TH>Variation</TH><TH>Min. Population</TH><TH>Max. Population</TH></TR>\n";
	ss<<"\t<TR><TD>"<<initPopulationSize<<"</TD><TD>"<<-1.0 * growthRate<<"</TD><TD>"<<carry<<"</TD>";
	ss<<"<TD>"<<timeOfMaxGrowth<<"</TD><TD>"<<polarity<<"</TD><TD>"<<variation<<"</TD><TD>"<<minPoolSize<<"</TD><TD>"<<maxPoolSize<<"</TD></TR></TABLE>\n";
	return ss.str();
}

string RichardsLogistic::GenerateCfgString() {
	stringstream ss;
	ss<<"GROWTH_RATE RICHARDS "<<initPopulationSize<<" "<<variation<<" "<<-1.0 * growthRate<<" "<<carry<<" "<<timeOfMaxGrowth<<" "<<polarity;
	return ss.str();
}


}

}
