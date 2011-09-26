//
// C++ Interface: population
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPOPULATION_H
#define SIMULATIONPOPULATION_H
#include "utility/types.h"
#include "chrompool.h"
#include "modelmanager.h"

namespace Simulation {

class PoolManager;

/**
 * @brief Base class for all sample types
 * Each Sample object represents a single group of datasets. They produce the datasets according to their configuration
 * and are capable of writing data in the correct format
 * Samples draw individuals from the pool manager- so more than 1 sample can be drawn from a set of chromosome pools
 * Status and is applied through the application of a model manager (which ultimately obtains a model). 
 */
class Sample {
public:
	Sample(float genoError, float phenoError, float missingData) : genoError(genoError), phenoError(phenoError), missingData(missingData) { }
	virtual ~Sample() {	}
	/**
	 * @brief Constructs a population based on the genetic data in pools with status derived from a 
	 * model provided by models.
	 * @param pools Array of chromosome pools. Each pool contains genetic data to be used to create a single individual
	 * @param chromCount Size of the array above
	 * @param models is used to randomly assign status
	 * @param individualCount Size of the resulting collection of individuals
	 */
	virtual void BuildSample(PoolManager &pools, ModelManager& models, uint individualCount)=0;

	/**
	 * @brief Adds individuals to the local sample based on the model passed
	 * @param pools Array of gene pools
	 * @param chromCount Size of the array above
	 * @param model the model used to apply status
	 * @param individualCount the number of individuals to be used
	 */
	virtual void BuildSample(PoolManager &pools, PenetranceModel* model, uint individualCount)=0;		

	/**
	 * @brief Apply genotype error based on each locus setting
	 */
	virtual void ApplyGenocopyError(PoolManager *mgr, uint familycount);
	
	/**
	 * @brief Apply phenocopy error 
	 */
	virtual void ApplyPhenocopyError(PoolManager *mgr, uint familycount)=0;

	/**
	 * @brief Apply a percentage to the dataset data and erase genotypes at the positions
	 */
	virtual void ApplyMissingData(PoolManager *mgr, uint familycount);

	/**
	 * @brief Works through the various types of errors and applies them each on the dataset
	 * according to their predefined settings
	 */
	virtual void ApplyErrors(PoolManager *mgr, uint familycount);
	/**
	 * @brief Resets the local data back to a starting point (for starting a new sample)
	 */
	virtual void Reset()=0;

	virtual void Purge();

	/**
	 * @brief Returns a simple segment used for generating file names. This will be something like: 15-85.cc for 15% affected case/control. The application is still required to generate the rest of the filename, however
	 */
	virtual string GetDescriptor()=0;

	/**
	 * @brief Dump the contents of the sample population to the stream
	 * @param os The stream to which the sample will be written
	 */
	virtual int WriteDataset(ostream &os, uint *gtCounts)=0;

	/**
	 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
	 * @param os The stream to which the sample will be written	
 	 */
	virtual void WritePhased(ostream &os)=0;

	virtual void GenerateReport(ostream &os, uint padding) = 0;

	virtual int GetSnpCount();
protected:
	vector<Individual*> people;				///<The actual population
	float genoError;						///<Genotype error associated with the local sample
	float phenoError;						///<Phenotype error
	float missingData;						///<Missing data percentage
};

/**
@brief Oversees the production of case/control datasets and saves them to file. 
These classes depend upon a remote GenePool array object to extract their individuals from

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BasicSample : public Sample {
public:
	/**
	 * @brief Construction
	 * @param percentAffected Percentage of individuals in the dataset to be found affected
	 * @param individualCount The total number of individuals to be found in a single dataset
	 */
    BasicSample(float percentAffected, float genoErr, float phenoErr, float missing);

    virtual ~BasicSample();

	/**
	 * @brief Constructs a population based on the genetic data in pools with status derived from a 
	 * model provided by models.
	 * @param pools Array of chromosome pools. Each pool contains genetic data to be used to create a single individual
	 * @param chromCount Size of the array above
	 * @param models is used to randomly assign status
	 * @param individualCount Size of the resulting collection of individuals
	 */
	virtual void BuildSample(PoolManager &pools, ModelManager& models, uint individualCount);

	/**
	 * @brief Adds individuals to the local sample based on the model passed
	 * @param pools Array of gene pools
	 * @param chromCount Size of the array above
	 * @param model the model used to apply status
	 * @param individualCount the number of individuals to be used
	 */
	virtual void BuildSample(PoolManager &pools, PenetranceModel* model, uint individualCount);		

	/**
	 * @brief Resets the local data back to a starting point (for starting a new sample)
	 */
	virtual void Reset();


	/**
	 * @brief Returns a simple segment used for generating file names. This will be something like: 15-85.cc for 15% affected case/control. The application is still required to generate the rest of the filename, however
	 */
	virtual string GetDescriptor();

	/**
	 * @brief Dump the contents of the sample population to the stream
	 * @param os The stream to which the sample will be written
	 */
	virtual int WriteDataset(ostream &os, uint *gtCounts);

	/**
	 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
	 * @param os The stream to which the sample will be written	
 	 */
	virtual void WritePhased(ostream &os);

	void GenerateReport(ostream &os, uint padding);

	void ApplyPhenocopyError(PoolManager *pools, uint familyCount);
	
protected:
	float percentAffected;					///<Percentage of overall individuals who are affected

};

inline
BasicSample::BasicSample(float percAffected, float genotypeError, float phenoError, float missing) : 
		Sample(genotypeError, phenoError, missing), percentAffected(percAffected) {
	assert(percAffected < 1.0 && percAffected > 0.0);
}

inline
BasicSample::~BasicSample() {
 }


inline
void BasicSample::Reset() {
	people.clear();
}

inline
int Sample::GetSnpCount() {
	int count = 0;
	
	if (people.size() > 0)
		count = people[0]->GetSnpCount();

	return count;
}

inline
string BasicSample::GetDescriptor() {
	stringstream ss;
	ss.precision(0);
	uint perc = (uint)(percentAffected *100.0);
	if (genoError > 0.0)
		ss<<"-"<<setprecision(2)<<(genoError*100.0)<<"GE";
	if (phenoError > 0.0)
		ss<<"-"<<setprecision(2)<<(phenoError*100.0)<<"PE";
	if (missingData > 0.0)
		ss<<"-"<<setprecision(2)<<(missingData*100.0)<<"M";
	ss<<setprecision(2)<<"-"<<perc<<"A"<<".mdr";
	return ss.str();
}

}

#endif
