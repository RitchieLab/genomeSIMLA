//
// C++ Implementation: populationcontroller
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "populationcontroller.h"

namespace GenomeSIM {

namespace GUI {

PopulationController::PopulationController() : 
			modelType(LinearGrowthRate), linearModel(10000, 0.0, 0.0), expModel(1500, 0.03, 0.0), 
			logModel(1500, 0.03, 120000, 0.0), 
			richModel(1500, 120000, 500, 0.02, 0.1, 0.0), 
			growthModel(&linearModel)	{
	rates = new GrowthRate*[4];
	rates[LinearGrowthRate] = &linearModel;
	rates[ExponentialGrowthRate] = &expModel;
	rates[LogisticGrowthRate] = &logModel;
	rates[RichardsLogisticGrowthRate] = &richModel;
}


PopulationController::~PopulationController()
{	

}

void PopulationController::SetCurrentGrowthModel(GrowthRate *model, GrowthRateType modelType) {
	this->modelType = modelType;

	if (modelType == LinearGrowthRate) 
		linearModel 	= *((LinearGrowth*)model);

	else if (modelType == ExponentialGrowthRate) 
		expModel 		= *((ExponentialGrowth*)model);

	else if (modelType == LogisticGrowthRate)
		logModel 		= *((LogisticGrowth*)model);

	else if (modelType == RichardsLogisticGrowthRate)
		richModel 		= *((RichardsLogistic*)model);

	RefreshControls(modelType);
}

string PopulationController::RefreshReport(const char *filename, int start, int stop, int interval) {
	assert(growthModel);
	return  growthModel->DrawGrowthChart( filename, start, stop, interval, 800, 400, true);
}


GrowthRate *PopulationController::GetGrowthModel() { 	
	return growthModel;
}


}

}
