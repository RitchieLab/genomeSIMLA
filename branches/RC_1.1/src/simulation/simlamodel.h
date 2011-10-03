//
// C++ Interface: simlamodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONSIMLAMODEL_H
#define SIMULATIONSIMLAMODEL_H

#include "penetrancemodel.h"
#include "individual.h"
#include "simla/penetrance.h"
#include "poolmanager.h"

namespace Simulation {
namespace StatusModel {

/**
Adaptor object for simla disease models

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SimlaModel : public PenetranceModel {
public:
    SimlaModel();

	/**
	 * @brief Initialize the Penetrance object for use.
	 * @param cfg This is the configuration file to be used by the simla component
	 * @param pools The pool manager used to generate the data for priming simla's algorithm
	 */
	//virtual void Initialize(char *cfg, PoolManager *pools);

	/**
	 * @brief The configuration of simla just tells how many and where the interactions occur. This specifies exactly which loci we are referring to
	 */
	//void AddDiseaseLoci(int chrom, int locus);

	void Load();

	bool Init(istream& i, PoolManager *pools);
	
	int BuildConvertedGenotypeIndex(vector<uint> & genotypes);

    ~SimlaModel();

	string Details();

	void GenerateReport(std::ostream &os, uint headerWidth);

	void GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci);

	bool ConfigureLocus(uint chrID, uint locID, float beta, float maf, float type, bool dxAtAl1);
	
	void SetupInteractions();
	void Refresh(PoolManager *pool);
	string GetModelConfiguration();
	//bool IsAffected(std::vector<uint>& genotypes);
	int GetType() { return 3; }
protected:
	int ConvertGenotype(int index, int genotype);
	/**
	 * @brief Converts multilocus genotypes indices into simla indices (required to re-order simla's table)
	 */
	int ConvertSimlaGenotype(int mlGenotype);
	struct ModelDetails {	
		float beta;
		float modelType;
		float dxFreq;
		string locusID;
		int idx;
		bool dxAtAllele1;

		ModelDetails() : 
				beta(0.0), modelType(0.0), dxFreq(0.0), idx(0), dxAtAllele1(false) { }
		ModelDetails(const ModelDetails& o) : 
				beta(o.beta), modelType(o.modelType), 
				dxFreq(o.dxFreq), locusID(o.locusID), idx(o.idx), dxAtAllele1(o.dxAtAllele1) { 
			cout<<"Simla Model Details: "<<beta<<" : "<<modelType<<" : "<<dxFreq<<" : "<<locusID<<" : "<<idx<<" : "<<dxAtAllele1<<"\n";
		}
		bool operator<(const ModelDetails& o) const { return beta<o.beta; } 
		ModelDetails &operator=(const ModelDetails o) { 
				beta=o.beta; modelType=o.modelType; dxFreq=o.dxFreq; 
				locusID=o.locusID; idx=o.idx; dxAtAllele1=o.dxAtAllele1; 
			cout<<"Simla Model Details: "<<beta<<" : "<<modelType<<" : "<<dxFreq<<" : "<<locusID<<" : "<<idx<<" : "<<dxAtAllele1<<"\n";
			return *this; }
		bool operator==(const ModelDetails& o) const { 
				return beta == o.beta && modelType == o.modelType && 
				dxFreq==o.dxFreq && dxAtAllele1==o.dxAtAllele1; } 
		ModelDetails(int idx, string locusID, float b, float m, float dxFreq, bool dxAtAllele1) : 
				beta(b), modelType(m), dxFreq(dxFreq), locusID(locusID), idx(idx), dxAtAllele1(dxAtAllele1) { 
			cout<<"Simla Model Details: "<<beta<<" : "<<modelType<<" : "<<dxFreq<<" : "<<locusID<<" : "<<idx<<" : "<<dxAtAllele1<<"\n";
		}
	};

	bool isInitialized;						///<Make sure we don't try to use this without priming the algorithm
	Simla::Penetrance penetranceModel;		///<This will work the magic that is simla
	float prevalence;						///<Prevalence
	map<string, ModelDetails> betaValues;	///<Beta values for the individual loci
	vector<string> betaKeys;				///<We need to know the order they were give to us
	int interactionCount;					///<The number of loci involved in interactions
	string configFilename;					///<Contains the details of the model
	int chromCount;
	PoolManager *poolMgr;

	void GenerateSeedData(int seedSize, int *seedData);
	void DrawIndividual(Individual &person);
	void ReserveIndividual(Individual &person);
	
};



inline
string SimlaModel::Details() {
	return string("Simla Model ");// + filename;
}

}

}

#endif
