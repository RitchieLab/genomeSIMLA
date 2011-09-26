//
// C++ Interface: diseasemodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONDISEASEMODEL_H
#define SIMULATIONDISEASEMODEL_H
#include "utility/types.h"

namespace Simulation {

/**
 * @brief Base class for disease models. 
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class DiseaseModel  {
public:
	/**
	 * @brief construction
	 * @param modelID The unique ID for the model	
 	 * @param prob The probability an individual will be affected with this type
	 */
	DiseaseModel(uint modelID, float prob) : modelID(modelID), probability(prob) {}
	virtual ~DiseaseModel() {}

	/**
	 * @brief Returns probability associated with the local model
	 */
	virtual float GetProbability() { return probability; } 

	/**
	 * @brief Returns the unique ID for the model 
	 */
	virtual uint GetModelID(uint modelID) { return modelID; }

	

	/**
	 * Performs a quick evaluation of the individual's genetic makeup and determines status randomly
	 * 
	 */
	virtual bool IsAffected(std::vector<uint>& genotypes) = 0;

	//Moved to samples so we could do datasets with these variations
	//static float phenocopyError;				///<Probability of affected status == true due to phenocopy error

protected:
	uint modelID;						///<Unique ID for a given model
	float probability;					///<This is used to determine percentage of a population to be assed by this model
};

}

#endif
