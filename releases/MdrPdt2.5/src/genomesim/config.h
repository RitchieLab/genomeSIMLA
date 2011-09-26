//
// C++ Interface: config
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONCONFIG_H
#define SIMULATIONCONFIG_H

#include "defaults.h"

#include "utility/utility.h"
#include "simulation/modelmanager.h"
#include "simulation/basicsample.h"
#include "simulation/poolmanager.h"
#include "simulation/growthrate.h"
#include "simulation/gamodel.h"

namespace GenomeSIM {

using namespace Simulation;

using namespace std;
using namespace Utility;
using namespace Simulation::PopulationGrowth;
/**
@brief The configuration parser for genomeSIM

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Config : public ConfigurationParser {
public:
    Config();

    ~Config();

	/**
 	 * Dumps a Text Report summarizing the configuration options encountered
	 */
	virtual void GenerateReport(ostream& os);	
	
	/**
	 * @brief Called to assign <i>value</i> to the entry, <i>key</i>
	 * @param key Left hand portion of the configuration line
	 * @param value Right hand portion of the configuration line
	 * @return true/false to indicate that it was a valid entry that was handled
	 */
	virtual bool SetValues(const char *key, const char *value, const char *remainder);

	/**
	 * @brief Check for missing data or impossible configuration settings
	 */
	virtual bool Validate();

	/**
	 * @brief This is called after everything is loaded to aggregate filenames and other composite information
	 */
	bool PostLoad();

	
	/**
	 * @brief This returns the application's pool manager already properly loaded
	 */
	PoolManager *GetPoolManager();

	/**
	 * @brief Returns the growth model described in the configuration file
	 */
	GrowthRate *GetGrowthRate();

	void SetMappingFn(const char *val);

	void LoadChromosome(const char *line);
	
	struct General {
		string outputName;					///<The first part of the reports. This is specified on the cmd line
		uint seed;							///<The random seed 
		uint poolSize;						///<The number interchangeable chromosomes in a given pool
		uint firstGeneration;				///<Which generation did we start with. Used when picking up where we left off
		uint firstDropPoint;				///<The last generation 
		uint dropFrequency;					///<how many generations between drops
		uint dropCount;						///<How many drops are expected be performed
		bool createLDMaps;					///<Do we want to create LD maps using haploview?
		bool createSampleLDmaps; 			///<Do we want to generate LD maps for our samples?
		bool keepPhasedOutput;				///<Do we want to keep the phased data after it is done?
		General() : outputName(""), seed(DEFAULT_SEED), poolSize(DEFAULT_POOL_SIZE), firstGeneration(0), firstDropPoint(1000), dropFrequency(250), dropCount(4), createLDMaps(false), createSampleLDmaps(false), keepPhasedOutput(true) {}
	} generalSettings;

	struct Status {
		//float phenoCopy;					///<fraction of affected that are due to environmental factors
		//float genotypeError;				///<The amount of error in genotypes
		float heteroRate;
		ModelManager models;
		bool doCreateModel;	
		void AddModel(const char *line);
		void DefineModel(const char *line, PoolManager *poolMgr);
		Status() : heteroRate(0.0) {}
	} statusSettings;

	struct Dataset {
		vector<Sample *> samples;			///<Describe the types of datasets to be used
		uint simSets;						///<The number of datasets per filetype to be written
		uint individualCount;				///<number of indidivuals to be written to a file (this might be interpreted differently for certain 
		bool writeDatasets;
		~Dataset() { 
			uint count = samples.size();
			for (uint i=0; i<count; i++)
				delete samples[i];
		}
		Dataset() : simSets(100), individualCount(100), writeDatasets(false) { }
	} datasetSettings;						///<Types of datasets

	struct ConfigKeywords {
		//General settings
		static char *RandomSeed;
		static char *PoolSize;
		static char *MaxPoolSize;
		static char *MinPoolSize;
		static char *AlleleLimits;
		static char *SetAlleleFrequency;
		static char *SetAlleleFreqRange;
		static char *GrowthRate;
		static char *DefaultAlleleFreq;
		
		//Status assignment
		static char *PhenoCopy;
		static char *GenotypeError;

		static char *FirstDropPoint;
		static char *DropFrequency;
		static char *DropCount;

		static char *AddModel;
		static char *DefineModel;
		static char *GenerateModels;
		static char *IndividualCount;
		static char *FamilyCount;			///<Number of families. This is actually the same as individual count
		static char *DatasetCount;	
		static char *StandardPedigreeHeader;
		
		static char *SetDefaultBlock;
		static char *AddChromosome;			///<Start the definition of a new chromosome
		static char *LoadChromosome;		///<Loads a definition from the chromosome
		static char *AddBlock;				///<Add a block to the current chromosome
		static char *RandomizeAlleleFreq;	///<Turns on/off randomization of allele frequencies
		static char *AddDataset;			
		static char *AddFamily;
		static char *AddFamilyRandom;
		static char *AddCC;

		static char *PathToJava;	
		static char *JavaSettings;
		static char *PathToHaploview;
		static char *HaploviewSettings;
		static char *ProduceLDMaps;
		static char *KeepPhasedOutput;
		static char *HaploviewOverviewSettings;
		static char *GenerateOverviewLD;
		static char *GenerateSampleLDMaps;
		static char *HaploviewWindowSize;
		static char *HaploviewWindowStride;
		static char *MappingFN;

		static char *GAFitnessThreshold;
		static char *GAIterations;

	};
	//string outputPrefix;					///<This is used in the generation of each filename
	string errorMsg;						///<Used to report any parsing errors;
protected:
	PoolManager pools;
	ChromPool::BlockDefinition defaultBlock;///<Used to define/store the default block qualities
	ChromPool *currPool;					///<Just used to keep up with the current pool. It should already be in the manager


	string GetRemainderString(const char *line);
	void SetAlleleFrequency(const char *line);
	void SetAlleleFreqRange(const char *line);
	void AddChromosome(const char *line);
	void AddBlock(const char *line);
	void SetDefaultBlock(const char *line);
	void SetGrowthRate(const char *line);
	void SetDefaultAlleleFreq(const char *line);
	void AddDataset(const char *line);

	GrowthRate *growthModel;
};

					

inline
void Config::Status::DefineModel(const char *line, PoolManager *poolMgr) {
	stringstream ss(line);
	string name;
	float percentage;
	string cmd;
	bool doContinue = true;
	int chrID, locID;					///<Used to determine disease loci
	
	vector<int> chrIDs;
	vector<int> locIDs;

	ss>>cmd>>name>>percentage;

	while (!ss.eof()) {
		ss>>chrID;
		if (!ss.eof()) {
			ss>>locID;

			chrIDs.push_back(chrID-1);
			locIDs.push_back(locID-1);
		}

		else {
			cout<<"Misconfigured model definition: \n\t"<<line<<"\n";
			cout<<"Correct Syntax: \tDEFINE_MODEL model_name probability chr1 locus1 [chrN locusN...]\n";
			abort();
		}
	}
	PenetranceModel *newModel = new GAModel(models.GetModelCount(), name.c_str(), percentage, chrIDs, locIDs);
	((GAModel*)newModel)->SetPoolManager(poolMgr);
	models.AddModel(newModel);
}



inline
void Config::Status::AddModel(const char *line) {
	stringstream ss(line);
	string filename = "";
	float percentage = 0.0;
	string cmd;

	ss>>cmd>>filename>>percentage;

	if (percentage > 0.0 ) { 
		PenetranceModel *newModel = new PenetranceModel(models.GetModelCount(), percentage);
		newModel->SetFilename(filename.c_str());
		//newModel->Load();
		models.AddModel( newModel);
	}
	else {
		cout<<"Misconfigured model definition: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \tADD_MODEL filename probability\n";
		abort();
	}
}

}

#endif
