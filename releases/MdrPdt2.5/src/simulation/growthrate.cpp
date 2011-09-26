//
// C++ Implementation: growthmodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "growthrate.h"
#include "exponentialgrowth.h"
#include "logisticgrowth.h"
#include "lineargrowth.h"

namespace Simulation {
namespace PopulationGrowth {

uint GrowthRate::maxPoolSize = 10000;		///<App wide max size. This is valid for ALL models
uint GrowthRate::minPoolSize = 10000;		///<App wide min
float GrowthRate::variation =  0.0005;	 	///<App wide variation to be applied to size calculation


void GrowthRate::GenerateReport( std::ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<"Fixed Population Size: " << maxPoolSize<<endl;
}

/*** No Growth *************************************************/
GrowthRate::GrowthRate() { }

GrowthRate::~GrowthRate() { }

uint GrowthRate::GetPopulationSize(uint t) { 
	return maxPoolSize;
}

uint GrowthRate::operator()(uint t) {
	return maxPoolSize;
}

double GrowthRate::AdjustedGrowth(double size) {
	//double adjustment = ((float)GetInitialPopulationSize() * rnd(variation)) - ((float)GetInitialPopulationSize()* variation/2.0);
	double adjustment = (size * rnd(variation)) - (size*variation/2.0);
	return size + adjustment;
}

uint GrowthRate::GetInitialPopulationSize() {
	return maxPoolSize;
}

string GrowthRate::GetType() {
	return "Fixed";
}

GrowthRate *GrowthRate::LoadModel(istream& ss) {
	GrowthRate *model = NULL;
	string type;
	ss>>type;
	
	if (strcmp(type.c_str(), "EXPONENTIAL") == 0) {
		float rate;
		uint initialPop;
		ss>>initialPop>>variation>>rate;
		if (rate > 0) 
			model = new ExponentialGrowth(initialPop, rate);
		else 
			model = new GrowthRate();

	}
	else if (strcmp(type.c_str(), "LOGISTIC") == 0) {
		uint carry;
		float rate; 
		uint initialPop;
		ss>>initialPop>>variation>>rate>>carry;
	
		if (rate > 0) 
			model = new LogisticGrowth(initialPop, rate, carry);
		else 
			model = new GrowthRate();
	}
	else if (strcmp(type.c_str(), "LINEAR") == 0) {
		float rate;
		uint initialPop;
		ss>>initialPop>>variation>>rate;

		if (rate > 0) 
			model = new LinearGrowth(initialPop, rate);
		else 
			model = new GrowthRate();
	}
	return model;
}

void GrowthRate::DiagramGrowth(std::ostream& os, int start, int stop, int interval) {
	os<<"Generation,Population Size\n";
	for (int i=start; i<=stop; i+=interval) 
		os<<i<<","<<GetPopulationSize(i)<<endl;
}


}

}
