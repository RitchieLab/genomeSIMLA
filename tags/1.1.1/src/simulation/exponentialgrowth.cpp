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
ExponentialGrowth::ExponentialGrowth(uint initSize, float rate, float var) : GrowthRate(var), initPopSize(initSize), rate(rate) {
	 modelType = ExponentialGrowthRate;
	r0=log(1+rate);
	cout<<"Expontial R0: "<<r0<<"\n";
}



string ExponentialGrowth::GetHtmlReport() {
	stringstream ss;
	ss<<"<TABLE><TR><TH>Initial Population</TH><TH>Growth Rate</TH><TH>Variation</TH><TH>Min. Population</TH><TH>Max. Population</TH></TR>\n";
	ss<<"\t<TR><TD>"<<initPopSize<<"</TD><TD>"<<rate<<"</TD><TD>"<<variation<<"</TD><TD>"<<minPoolSize<<"</TD><TD>"<<maxPoolSize<<"</TD></TR></TABLE>\n";
	return ss.str();
}

uint ExponentialGrowth::GetPopulationSize( uint generation ) {
	double basePop = initPopSize * exp(r0 * (double)generation);
	double newSize = AdjustedGrowth(basePop);

	if (newSize == HUGE_VAL)
		newSize = maxPoolSize;
	if (newSize > maxPoolSize)
		newSize = maxPoolSize;
	else if (newSize < minPoolSize)
		newSize = minPoolSize;


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
	os<<setw(headerWidth)<<"Variation: "<<setprecision(4)<<variation<<endl;}

string ExponentialGrowth::GenerateCfgString() {
	stringstream ss;
	ss<<"GROWTH_RATE EXPONENTIAL "<<initPopSize<<" "<<variation<<" "<<rate;
	return ss.str();
}
}

}
