/* 
 * File:   basicheterogeneity.h
 * Author: torstees
 *
 * Created on January 13, 2011, 12:44 PM
 */

#ifndef BASICHETEROGENEITY_H
#define	BASICHETEROGENEITY_H

#include "penetrancemodel.h"
#include <vector>
#include <string>
#include <set>

namespace Simulation {
namespace StatusModel {

class BasicHeterogeneity : public PenetranceModel {
public:
	BasicHeterogeneity();
	BasicHeterogeneity(const BasicHeterogeneity& orig);
	virtual ~BasicHeterogeneity();

	/**
	 * @brief Load the model at the local file, filename
	 */
	virtual void Load();

	void Refresh(PoolManager *pools);

	void AddModel(const char *cfg);

	virtual void GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci);

	virtual bool Init(istream& s, PoolManager *pools);
	//bool Magic(PoolManager* pools);
	virtual string GetModelConfiguration();

	static PenetranceModel *GetModel(istream& details, PoolManager *pools);

	int GetStatus(std::vector<uint>& genotypes, float& outcome);

	std::vector<uint> GetLocalGenotypes(PenetranceModel* model);

	void GenerateReport(std::ostream &os, uint headerWidth);

	std::string Details();
private:
	std::vector<std::string> modelConfigurations;
	std::set<DiseaseLocus> allLoci;
	std::vector<PenetranceModel *> models;
	std::map<PenetranceModel *, std::vector<int> > modelLoci;
	int totalCaseCount;
};

#endif	/* BASICHETEROGENEITY_H */





}

}
