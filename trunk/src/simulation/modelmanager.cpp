//
// C++ Implementation: modelmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <assert.h>
#include "utility/random.h"
#include "modelmanager.h"

namespace Simulation {

using namespace std;

ModelManager::ModelManager() : totalFrequency(0.0) {
}


ModelManager::~ModelManager() {
	Reset();
}

uint ModelManager::GetModelCount() {
	return models.size();
}

PenetranceModel *ModelManager::GetModel(uint idx) {
	assert(idx<models.size());
	return models[idx];
}

PenetranceModel *ModelManager::GetModel() {
	double roll = Utility::Random::globalGenerator.drand();
	double currValue = 0.0;
	PenetranceModel *selectedModel = NULL;
	uint modelCount = models.size();
	
	for (uint i=0; i<modelCount && !selectedModel; i++) {
		currValue += models[i]->GetProbability();
		if (roll <= currValue) 
			selectedModel=models[i];
	}
	return selectedModel;	
}

bool ModelManager::Validate() {
	return totalFrequency > 0.999 && totalFrequency < 1.001;
}
		

void ModelManager::GenerateReport(ostream &os, uint width) {
	uint count = GetModelCount();
	for (uint i=0; i<count; i++) {
		models[i]->GenerateReport(os, width);
		os<<endl;
	}
}

void ModelManager::LoadModelData() {
	uint count = GetModelCount();
	for (uint i=0; i<count; i++) {
		models[i]->Load();
	}
}

void ModelManager::Reset() {
	totalFrequency = 0.0;
	vector<PenetranceModel *>::iterator i = models.begin();
	vector<PenetranceModel *>::iterator end = models.end();
	
	for (; i!=end; i++) 
		delete (*i);
	models.clear();
}

bool ModelManager::AddModel( PenetranceModel *model) {
	bool success = false;
	if (totalFrequency + model->GetProbability() <= 1.0) {	
		totalFrequency+=model->GetProbability();
		models.push_back(model);
		success=true;
	}
	return success;
}

}
