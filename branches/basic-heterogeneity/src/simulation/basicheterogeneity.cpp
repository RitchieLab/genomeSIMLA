/* 
 * File:   basicheterogeneity.cpp
 * Author: torstees
 * 
 * Created on January 13, 2011, 12:44 PM
 */

#include "basicheterogeneity.h"

namespace Simulation {
namespace StatusModel {


BasicHeterogeneity::BasicHeterogeneity() {
}

BasicHeterogeneity::BasicHeterogeneity(const BasicHeterogeneity& orig) {
}

BasicHeterogeneity::~BasicHeterogeneity() {
}

void BasicHeterogeneity::Load() {
	std::vector<PenetranceModel*>::iterator itr = models.begin();
	std::vector<PenetranceModel*>::iterator end = models.end();

	while (itr != end) {
		(*itr++)->Load();
	}
}

void BasicHeterogeneity::Refresh(PoolManager *pools) {
	size_t lociCount = loci.size();

	for (size_t i=0; i<lociCount; i++) {
		DiseaseLocus &locus = loci[i];
		pools->GetAlleleFrequency(locus.chromosome, locus.locusIdx, locus.alFreq1, locus.alFreq2);
	}
}
std::string BasicHeterogeneity::GetModelConfiguration() {
	return "This hasn't been set up for Basic Heterogeneity yet....sorry\n";
}

bool BasicHeterogeneity::Init(istream& s, PoolManager *pools) {
	bool success = true;
	std::vector<std::string>::iterator itr = modelConfigurations.begin();
	std::vector<std::string>::iterator end = modelConfigurations.end();
	models.clear();

	while (itr != end) {
		stringstream ss(itr++->c_str());
		string cmd;

		ss>>cmd;
		Simulation::StatusModel::PenetranceModel *model = Simulation::StatusModel::PenetranceModel::GetModel(ss, pools, modelConfigurations);
		if (model == NULL)
			return false;
		
		int modelSize = model->GetModelSize();
		for (int i=0; i<modelSize; i++) {
			DiseaseLocus locus = ((DiseaseModel*)model)->GetLocus(i);
			allLoci.insert(locus);
		}
		models.push_back(model);
		modelLoci[model] = vector<int>();
	}

	SetModelSize(allLoci.size());

	std::set<DiseaseLocus>::iterator lItr = allLoci.begin();
	std::set<DiseaseLocus>::iterator lEnd = allLoci.end();

	while (lItr != lEnd) {
		DiseaseLocus l = *lItr;
		AddDiseaseLoci(l.label.c_str(), l.chromosome, l.locusIdx, l.alFreq1);
		lItr++;
	}

	std::map<PenetranceModel*, vector<int> >::iterator mItr = modelLoci.begin();
	std::map<PenetranceModel*, vector<int> >::iterator mEnd = modelLoci.end();
	int locusCount = loci.size();
	while (mEnd != mItr) {
		vector<DiseaseLocus> localModel = mItr->first->GetDiseaseLoci();
		for (int i=0; i<locusCount; i++) {
			if (std::find(localModel.begin(), localModel.end(), loci[i])!=localModel.end())
				mItr->second.push_back(i);
		}

		mItr++;
	}

	return true;
}

std::string BasicHeterogeneity::Details() {
	return "Heterogenity model";
}
void BasicHeterogeneity::GenerateReport(std::ostream &os, uint headerWidth) {
	os<<"Heterogeneity Model\n";
}
void BasicHeterogeneity::GenerateDetailedReport(ostream& os, vector<Locus*> &diseaseLoci) {
	std::vector<PenetranceModel*>::iterator itr = models.begin();
	std::vector<PenetranceModel*>::iterator end = models.end();

	while (itr != end) {
		(*itr++)->GenerateDetailedReport(os, diseaseLoci);
	}
}

void BasicHeterogeneity::AddModel(const char *cfg) {
	modelConfigurations.push_back(cfg);
}

PenetranceModel *BasicHeterogeneity::GetModel(istream& ss, PoolManager *pools) {
	PenetranceModel *model = new BasicHeterogeneity();
	return model;
}


int BasicHeterogeneity::GetStatus(std::vector<uint>& genotypes, float& outcome) {
	PenetranceModel *model = models[totalCaseCount % models.size()];
	std::vector<uint> localGenotypes;

	std::vector<int> &indices = modelLoci[model];
	std::vector<int>::iterator itr = indices.begin();
	std::vector<int>::iterator end = indices.end();

	while (itr != end) {
		localGenotypes.push_back(genotypes[*itr++]);
	}
	uint status = model->GetStatus(localGenotypes, outcome);
	if (status) {
		totalCaseCount++;
	}
	return status;
}

}
}

