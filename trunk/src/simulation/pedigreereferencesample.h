//
// C++ Interface: pedigreereferencesample
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_PEDIGREETEMPLATESPEDIGREEREFERENCESAMPLE_H
#define SIMULATION_PEDIGREETEMPLATESPEDIGREEREFERENCESAMPLE_H

#include "templatedpedigree.h"
#include "pedigreesample.h"
namespace Simulation {


/**
 * Generates a template and upon which the pedigree structures of the datasets will be based.
	@author Eric Torstenson 
*/
class PedigreeReferenceSample : public PedigreeSample {
public:
    PedigreeReferenceSample(int replicates, float gtError, float phError, float missingData, const char *desc);
    ~PedigreeReferenceSample();

	/**
	 * @brief Apply disease to the pedigree
	 */
	void BuildSample(PoolManager &pools, PenetranceModel* model);
	
	void ApplyPhenocopyError(PoolManager *pools, PenetranceModel *model);
	
	string GetDescriptor();
	
	string GetType() { return "PED_REFERENCE"; }
	
	string GetSummary();
	
	string GetDetails();
		
	string GetLabel() { return description; }	///<Quick description of the dataset
	
	string Details();							///<Return model details in string format
	
	bool UseMeta() { return true; }				///<Indicate that we actually can use meta data

	/**
	 * @brief Determine the number of replicates to be written to file
	 */
	void SetReplicateCount(int replicateCount);

	int WriteDataset(const char *baseFilename, uint *gtCount);

	void PrepSample(PoolManager *pools, PenetranceModel* model);
	void WriteLOD(const char *baseFilename);
	void SetReference(const char *filename);
	int WriteDataset(ostream &os, uint *gtCount);
	static float minMaxTolerance;	
	static uint maxAttemptsAtReplication;
	static bool DoWriteLOD;
	uint GetCount();							///< Returns the number of individuals in the reference set
protected:
	/**
	 * @brief these are the pedigree templates that represent our target disease model
	 */
	PedigreeTemplates::TemplatedPedigree pedigrees;
	int replicates;								///< The number of same pedigrees with varied status applied
	StatusModel::PenetranceModel *curModel;		///< We are capturing the disease model for replicating
	string reference;
	PedigreeTemplates::TemplatedPedigree::StatusCounts Replicate();
	PoolManager *pools;
	uint minAffectedCount;
	uint maxAffectedCount;
	bool writeGenotypesForAll;					///< Allow the user to override the normal behavior of only writing genotypes for persons of interest
};




}

#endif
