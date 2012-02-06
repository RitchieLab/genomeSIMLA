//
// C++ Interface: modelmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONMODELMANAGER_H
#define SIMULATIONMODELMANAGER_H

#include <iostream>
#include <vector>
#include "penetrancemodel.h"

namespace Simulation {

using namespace StatusModel;

/**
@brief Helps keep up with the various models and their frequencies

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ModelManager{
public:
    ModelManager();

    ~ModelManager();

	/**
	 * @brief Add a model to the manager. This assumes responsibility for the memory
	 */
	bool AddModel(PenetranceModel *model);

	/**
	 * @brief Returns a random model based on frequencies
	 */	
	PenetranceModel *GetModel();

	/**
	 * @brief returns the number of models being managed
	 */
	uint GetModelCount();
	
	/**
	 * @brief Get a specific model	
	 */
	PenetranceModel *GetModel(uint idx);

	/**
	 * @brief Either total frequencies should be near 1.0 or 0.0
	 */
	bool Validate();

	/**
	 * @brief clear all models and recent the total frequency
	 */
	void Reset();

	/**
	 * @brief Send configuration information to the stream
	 * @param os The stream to be written to
	 * @param width the width of the header column
	 */
	void GenerateReport(std::ostream &os, uint width);

	/**
	 * @brief Instructs each of the local models to load their contents
	 */
	void LoadModelData();
protected:
	float totalFrequency;						///<Frequency assocaited with all of the models
	std::vector<PenetranceModel *> models;		///<The collection of model objects
	
};


}

#endif
