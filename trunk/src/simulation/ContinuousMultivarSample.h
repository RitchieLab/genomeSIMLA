/*
 * ContinuousMultivarSample.h
 *
 *  Created on: Sep 7, 2012
 *      Author: jrw32
 */

#ifndef SIMULATION_CONTINUOUSMULTIVARSAMPLE_H
#define SIMULATION_CONTINUOUSMULTIVARSAMPLE_H

#include <string>
#include <map>
#include <set>
#include "basicsample.h"

#include <gsl/gsl_matrix.h>

using std::string;
using std::map;
using std::set;
using std::pair;

namespace Simulation {

class ContinuousMultivarSample : public Sample {
public:
	ContinuousMultivarSample(float genoError, float phenoError, float missingData, int n_in, const string& e_fn, const string& c_fn, const string& name_in);
	virtual ~ContinuousMultivarSample();

	virtual void BuildSample(PoolManager& pools, PenetranceModel* model);
	virtual void ApplyPhenocopyError(PoolManager* mgr, PenetranceModel* model);
	virtual void Reset();
	virtual void Purge();

	virtual bool requireModel() { return false; }

	virtual string Details();
	virtual string GetDescriptor();
	virtual string GetType();
	virtual string GetLabel();
	virtual string GetSummary();
	virtual string GetDetails();

	virtual int WriteDataset(std::ostream&, uint*);
	virtual int WriteBinaryDataset(std::ostream&, std::ostream&);
	virtual void GenerateReport(std::ostream&, uint);

	virtual void LoadBinarySample(ifstream* genotypes, ifstream* meta, uint peopleCount, uint chromCount, vector<uint>& chrId);
	virtual void WriteSnpDetails(const char*, PoolManager* pools);

private:
	// No copying or assignment (I hope)
	ContinuousMultivarSample(const ContinuousMultivarSample& other);
	ContinuousMultivarSample operator=(const ContinuousMultivarSample& other);

	gsl_matrix* covar_matrix;

	// Tne number of effects to simulate
	int n_effects;

	// This is a sparse representation of a matrix indexed by both effect and SNP

	// Because SNPs can be referred to by both their index AND their label, we
	// need to store both possible representations, but we can't resolve them
	// until we actually build the sample.  grrr....
	map<int, map<pair<int,int>, double> > effect_matrix_idx;
	map<int, map<string, double> > effect_matrix_lbl;

	const string effect_fn;
	const string covar_fn;
	const string name;

	// The number of people in the sample
	int n_indiv;

	// The actual effects we are simulating
	map<Individual*, vector<double> > effect_map;




};

} /* namespace Simulation */
#endif /* SIMULATION_CONTINUOUSMULTIVARSAMPLE_H */
