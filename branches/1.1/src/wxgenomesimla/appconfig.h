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
#include "simulation/simlamodel.h"
#include "simulation/locselection.h"

namespace GenomeSIM {

using namespace Simulation;

using namespace std;
using namespace Utility;
using namespace Simulation::PopulationGrowth;


typedef map<string, Simulation::Visualization::LocusSelection> LocusSelectionMap;

/**
@brief The configuration parser for genomeSIM

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class AppConfig : public ConfigurationParser {
public:
    AppConfig();

    ~AppConfig();

	/**
	 * @lets give it a pseudo singleton ability.
	 */
	static AppConfig *GetPrimaryInstance() { 
		return primaryInstance;
	}
	static void SetPrimaryInstance(AppConfig *cfg) {
		primaryInstance = cfg;
	}

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

	/**
	 * @brief Setup the location mapping function (kosambi/haldane)
	 */
	void SetMappingFn(const char *val);

	/**
	 * @brief Loads the chromosome (loci positions)
	 */
	void LoadChromosome(const char *line);
	
	/**
	 * @brief Loads Initiates the parsing of the configuration file
	 */
	bool LoadSettings(const char *filename);

	/**
	 * @brief Random seed
	 */
	void SetRandomSeed(uint seed);

	/**
	 * @brief returns the filename associated with the growth rate
	 */
	string GetGrowthChartFilename() { return generalSettings.growthChartFilename; }

	bool WriteConfiguration(const char *config);
	
	struct General {
		string outputName;					///<The first part of the reports. This is specified on the cmd line
		uint seed;							///<The random seed 
		uint poolSize;						///<The number interchangeable chromosomes in a given pool
		uint firstGeneration;				///<Which generation did we start with. Used when picking up where we left off
		uint firstDropPoint;				///<The last generation 
		uint dropFrequency;					///<how many generations between drops
		uint dropCount;						///<How many drops are expected be performed
		bool createLDMaps;					///<Do we want to create LD maps using haploview?
		bool writeDetailedLD;				///<Do we want to create Detailed (complete) LD plots?
		bool createSampleLDmaps; 			///<Do we want to generate LD maps for our samples?
		bool keepPhasedOutput;				///<Do we want to keep the phased data after it is done?
		bool doLoad;						///<Are we loading a previous generation? 
		bool closePoolsBetweenAdvancing;	///<Do we keep the pools open between advancements?
		string growthChartFilename;			///<The name of the file containing growth details

		/**
		 * @brief build up the selectors so we can have them ready when the chromosome pools are ready
	 	 */
		LocusSelectionMap locusSelections;
		vector<string> selectionRegions;	///<cache the selections until we can translate them into loci

		General();
		
		/**
	 	 * @brief Add a new locus selection object
	 	 */
		void AddSelector(const char *selLine);
		/**
	 	 * @brief Add a region to search within
	 	 */
		void AddSelectorRegion(const char *selRegion);


		void WriteLocusSelections(ostream &os);

		/**
		 * @brief 
		 */
		void ReconcileSelections(PoolManager &pools);
		void PartiallyReconcileSelections(PoolManager &pools);
	} generalSettings;

	struct Status {
		//float phenoCopy;					///<fraction of affected that are due to environmental factors
		//float genotypeError;				///<The amount of error in genotypes
		float heteroRate;
		//ModelManager models;
		PenetranceModel *model;
		bool doCreateModel;	
		string modelCfg;
		void SetModel(const char *line);
		void WriteModelSelection(ostream &file);

		PenetranceModel *LoadModel(const char *project, PoolManager *poolMgr);
		PenetranceModel *InitModel(PoolManager *poolmgr, bool requireModel);
		//void AddModel(const char *line);
		void DefineModel(const char *line);
		Status() : heteroRate(0.0), model(NULL) {}
		~Status() { delete model; }
	} statusSettings;

	struct Dataset {
		vector<Sample *> samples;			///<Describe the types of datasets to be used
		uint simSets;						///<The number of datasets per filetype to be written
		uint individualCount;				///<number of indidivuals to be written to a file (this might be interpreted differently for certain 
		bool writeDatasets;
		bool binaryDatasets;
		void Purge() {
			uint count = samples.size();
			for (uint i=0; i<count; i++)
				delete samples[i];
			samples.clear();
		}
		string GetDatasetConfiguration();
		void WriteDatasetSelections(ostream &file);
		~Dataset() { 
			Purge();
		}
		Dataset() : simSets(100), individualCount(100), writeDatasets(false), binaryDatasets(false) { }
	} datasetSettings;						///<Types of datasets

	struct AppConfigKeywords {
		//General settings
		static const char *RandomSeed;
		static const char *PoolSize;
		static const char *MaxPoolSize;
		static const char *MinPoolSize;
		static const char *AlleleLimits;
		static const char *SetAlleleFrequency;
		static const char *SetAlleleFreqRange;
		static const char *GrowthRate;
		static const char *DefaultAlleleFreq;
		
		//Status assignment
		static const char *PhenoCopy;
		static const char *GenotypeError;

		static const char *FirstDropPoint;
		static const char *DropFrequency;
		static const char *DropCount;

		//static char *AddModel;
		static const char *SetModel;
		static const char *DefineModel;
		static const char *LabelBased;
		static const char *IndexBased;
		static const char *Simla;
		static const char *Simpen;
		static const char *PenTable;
		static const char *GenerateModels;
		static const char *IndividualCount;
		static const char *FamilyCount;			///<Number of families. This is actually the same as individual count
		static const char *DatasetCount;	
		static const char *StandardPedigreeHeader;
		
		static const char *SetDefaultBlock;
		static const char *AddChromosome;			///<Start the definition of a new chromosome
		static const char *LoadChromosome;		///<Loads a definition from the chromosome
		static const char *SeedChromosome;		///<Creates a chromosome pool based on the seed data
		static const char *AddBlock;				///<Add a block to the current chromosome
		static const char *RandomizeAlleleFreq;	///<Turns on/off randomization of allele frequencies
		static const char *AddDataset;			
		static const char *AddFamily;
		static const char *AddFamilyRandom;
		static const char *AddFamilyType;
		static const char *AddCC;
		static const char *OffspringPerMating;	///<Keyword to set min/max children per mating (in advancement)

		static const char *GenerateOverviewLD;
		static const char *MappingFN;
		static const char *BlockReportSize;
		static const char *MaxSnpsPerRow;
		static const char *MaxSNPDistance;		///<For LD Comparisons

		static const char *GAFitnessThreshold;
		static const char *GAIterations;

		static const char *MinAlleleFreqThresh;
		static const char *IncludeAppConfiguration;

		static const char *DoDumpGenerationZero;

		static const char *TrueTypeFontName;		///<The font to be used
		static const char *BlockReportBuffer;		///<The number of snps to be used on either side of a block in an LD plot
		static const char *CssFilename;			///<Filename to be used as the CSS style sheet. This will be written to all HTML report pages and must be available in the location specified relative to the report whenever it is viewed
		static const char *ClosePoolsBetweenAdvancements;

		static const char *UseOriginalCrossing;
		static const char *PedUseOriginalCross;
		static const char *PhasedPedigrees;

		static const char *UseAltLD;
		static const char *DrawRSquared;
		static const char *DrawDPrime;
		static const char *WriteLDReport;
		//These two parameters allow for limiting the LD calculations in order to perform
		//scans in reasonable amounts of time
		static const char *MaxLdIndividualCount;	///<Number of individuals 
		static const char *PlotScanSize;			///<Number of snps

		static const char *BinaryDatasets;

		static const char *LocSelection;			///<Add a locus selection object
		static const char *AddRegion;				///<Add a region to a selection object

		static const char *SimultaneousChroms;	///<Number of chroms to be worked on at once
		static const char *ThreadsPerChrom;		///<Number of threads to be launched per chrom

		static const char *MaxThreadCount;		///<Specify the number of threads to be used
		static const char *OverideRecThreadCount;

		static const char *NoFastLD;				///<Turn off sampled LD
		static const char *TargetPopSize;			///<Pop. Size where GS will stop once it is reached
		static const char *LociPerChromReported;	///<Max Number of loci / chromosomes to be reported

		
	};
	//string outputPrefix;					///<This is used in the generation of each filename
	string errorMsg;						///<Used to report any parsing errors;

	string GetAppConfigurationFilename() { return configurationFilename; }
	void SetAppConfigurationFilename(const char *filename) { configurationFilename = filename; }

	void SetGrowthRate(const char *line);

	uint GetRandomSeed();

	ChromPool::BlockDefinition *GetDefaultBlock() { return &defaultBlock; }
	void ClearDatasets();
	void AddDataset(const char *line);
	void SetDefaultBlock(ChromPool::BlockDefinition& block);

	void ConfigureDiseaseModel(bool writeDatasets);

	void SummarizeDiseaseModel(ostream &os, vector<Locus *>diseaseLoci);
protected:
	PoolManager pools;
	ChromPool::BlockDefinition defaultBlock;///<Used to define/store the default block qualities
	ChromPool *currPool;					///<Just used to keep up with the current pool. It should already be in the manager


	string GetRemainderString(const char *line);
	void SetAlleleFrequency(const char *line);
	void SetAlleleFreqRange(const char *line);
	void AddChromosome(const char *line);
	void SeedChromosome(const char *line);
	void AddBlock(const char *line);
	void SetDefaultBlock(const char *line);
	void SetDefaultAlleleFreq(const char *line);
	
	bool doDumpGenerationZero;

	
	GrowthRate *growthModel;
	bool growthModelAppConfigured;
	static AppConfig *primaryInstance;

	string configurationFilename;
};

					

inline
void AppConfig::Status::DefineModel(const char *line) {
	modelCfg = line;
}

inline
void AppConfig::General::AddSelectorRegion( const char *selLine) {
	selectionRegions.push_back(selLine);
}

inline
AppConfig::General::General() : outputName(""), seed(DEFAULT_SEED), poolSize(DEFAULT_POOL_SIZE), 
		firstGeneration(0), firstDropPoint(1000), dropFrequency(250), dropCount(4), 
		createLDMaps(false), writeDetailedLD(false), createSampleLDmaps(false), 
		keepPhasedOutput(true), doLoad(false), closePoolsBetweenAdvancing(false) {

	//This default selection is just a place holder for "global" regions
	locusSelections["ALL_SELECTORS"] = Simulation::Visualization::LocusSelection(0.5, 0.0, 0.5, 10, 0, 50);
	locusSelections["ALL_SELECTORS"].SetLabel("ALL_SELECTORS");
}


/*
inline
void AppConfig::Status::SetModel(const char *line) {
	stringstream ss(line);
	string filename = "";
	float percentage = 0.0;
	string cmd;
	
	ss>>cmd>>filename>>percentage;

	if (model)  {
		cout<<"Overwriting previous model definition. (";
		model->GenerateReport(cout, 0);
		cout<<")\n";
		delete model;
	}

	if (percentage > 0.0) {
		model = new PenetranceModel();
		model->SetFilename(filename.c_str());
	}
}*/

/*
inline
void AppConfig::Status::AddModel(const char *line) {
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
*/
}

#endif
