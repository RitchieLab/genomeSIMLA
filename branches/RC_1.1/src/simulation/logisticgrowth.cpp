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
	os<<setw(headerWidth)<<"Carrying Capacity: "<<(uint)carry<<endl;
	os<<setw(headerWidth)<<"Variation: "<<setprecision(4)<<variation<<endl;
}




string LogisticGrowth::GetHtmlReport() {
	stringstream ss;
	ss<<"<TABLE><TR><TH>Initial Population</TH><TH>Growth Rate</TH><TH>Carrying Capacity</TH><TH>Variation</TH><TH>Min. Population</TH><TH>Max. Population</TH></TR>\n";
	ss<<"\t<TR><TD>"<<initPopulationSize<<"</TD><TD>"<<growthRate
		<<"</TD><TD>"<<carry<<"</TD><TD>"<<variation<<"</TD><TD>"<<minPoolSize<<"</TD><TD>"<<maxPoolSize<<"</TD></TR></TABLE>\n";
	return ss.str();
}




uint LogisticGrowth::FindInflection(uint startingPop, int  maxTries) {
	int x1=startingPop, x2 = 0;
	int y1=0, y2=0;
	y1 = (int)(carry/(1.0 + (C * exp (-1.0 * x1 * growthRate))));
	double lastGrowthRate = 0.0;
	bool doContinue = true;
	int n = 0;
	while (doContinue && n++ < maxTries) {
		y2 = (int) (carry/(1.0 + (C * exp (-1.0 * ++x2 * growthRate))));
		double growthRate = (double)(y2 - y1) / (double)(x2 - x1);
		doContinue = growthRate > 0.000001;
		//cout<<n<<"\t"<<y2<<"\t"<<growthRate<<"\n";
		//doContinue = lastGrowthRate < 0.0001 || lastGrowthRate > growthRate;
		lastGrowthRate = growthRate;
		x1=x2;
		y1=y2;
	}
	return x1;
}

string LogisticGrowth::GenerateCfgString() {
	stringstream ss;
	ss<<"GROWTH_RATE LOGISTIC "<<initPopulationSize<<" "<<variation<<" "<<growthRate<<" "<<carry;
	return ss.str();
}


}

}
