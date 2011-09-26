//
// C++ Interface: populationcontroller
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIPOPULATIONCONTROLLER_H
#define GENOMESIM_GUIPOPULATIONCONTROLLER_H

#include <string>
#include "utility/types.h"
#include "simulation/exponentialgrowth.h"
#include "simulation/lineargrowth.h"
#include "simulation/richardslogistic.h"

namespace GenomeSIM {

namespace GUI {

using namespace std;
using namespace Simulation::PopulationGrowth;



/**
&brief Responsible for all population growth interactions

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PopulationController{
public:


    PopulationController();

    virtual ~PopulationController();

	GrowthRate *GetGrowthModel();
	virtual string GetGrowthConfiguration() = 0;

	string RefreshReport(const char *filename, int start, int stop, int interval);

	void SetCurrentGrowthModel(GrowthRate *model, GrowthRateType modelType);

	virtual void Init() = 0;
protected:
	virtual void RefreshControls(GrowthRateType control) = 0;

	string growthConfiguration;
	GrowthRateType modelType;

	LinearGrowth linearModel;
	ExponentialGrowth expModel;
	LogisticGrowth logModel;
	RichardsLogistic richModel;
	GrowthRate *growthModel;

	GrowthRate **rates;							///<This is an index into each of the rates

};

}

}

#endif
