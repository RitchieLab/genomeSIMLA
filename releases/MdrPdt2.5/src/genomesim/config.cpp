//
// C++ Implementation: config
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <iostream>
#include <iomanip>
#include "config.h"
#include "utility/utility.h"
#include <assert.h>
#include <sstream>
#include "simulation/diseasemodel.h"
#include "simulation/pedigreesample.h"
#include "simulation/gamodel.h"

namespace GenomeSIM {

using namespace std;

char *Config::ConfigKeywords::PathToJava		= "PATH_TO_JAVA";
char *Config::ConfigKeywords::PathToHaploview	= "PATH_TO_HAPLOVIEW";
char *Config::ConfigKeywords::HaploviewSettings	= "HAPLOVIEW_SETTINGS";
char *Config::ConfigKeywords::ProduceLDMaps		= "PRODUCE_LD_MAPS";
char *Config::ConfigKeywords::KeepPhasedOutput	= "KEEP_PHASED";
char *Config::ConfigKeywords::HaploviewOverviewSettings = "HAPLOVIEW_OVERVIEW_SETTINGS";
char *Config::ConfigKeywords::GenerateOverviewLD = "GENERATE_OVERVIEW_LD";
char *Config::ConfigKeywords::HaploviewWindowSize = "HAPLOVIEW_WINDOW_SIZE";
char *Config::ConfigKeywords::HaploviewWindowStride = "HAPLVIEW_WINDOW_STRIDE";
char *Config::ConfigKeywords::GAFitnessThreshold	= "GA_FITNESS_THRESHOLD";
char *Config::ConfigKeywords::GAIterations 			= "GA_ITERATIONS";
char *Config::ConfigKeywords::GenerateSampleLDMaps = "GENERATE_SAMPLE_LD";
char *Config::ConfigKeywords::JavaSettings			= "JAVA_SETTINGS";
char *Config::ConfigKeywords::MappingFN			= "SET_MAPPING_FN";

////General Settings
char *Config::ConfigKeywords::PoolSize			= "POOL_SIZE";					///<Number of unique expressions of each of the blocks to be found in the various pools
char *Config::ConfigKeywords::MaxPoolSize		= "MAX_POOL_SIZE";				///<Ceiling for growth
char *Config::ConfigKeywords::MinPoolSize		= "MIN_POOL_SIZE";				///<Floor
char *Config::ConfigKeywords::RandomSeed		= "SEED";						///<Random seed
char *Config::ConfigKeywords::AlleleLimits		= "ALLELE_LIMITS";				///<Range of allele limits
char *Config::ConfigKeywords::GrowthRate		= "GROWTH_RATE";				///<Setup limits for growth rate / generation
char *Config::ConfigKeywords::DefaultAlleleFreq = "DEFAULT_ALLELE_FREQ";		///<The default frequencies of alleles
char *Config::ConfigKeywords::SetAlleleFrequency = "ALLELE_FREQUENCY";			///<Force a frequency on a given locus
char *Config::ConfigKeywords::SetAlleleFreqRange = "ALLELE_FREQ_RANGE";

////Status Settings
char *Config::ConfigKeywords::GenotypeError		= "GENOTYPE_ERROR";
char *Config::ConfigKeywords::PhenoCopy 		= "PHENOCOPY_ERROR";
char *Config::ConfigKeywords::AddModel			= "ADD_MODEL";
char *Config::ConfigKeywords::DefineModel		= "DEFINE_MODEL";
char *Config::ConfigKeywords::GenerateModels	= "GENERATE_MODELS";

/////Dataset settings
char *Config::ConfigKeywords::IndividualCount	= "INDIVIDUAL_COUNT";			///<Number of individuals per file (for Case/control)
char *Config::ConfigKeywords::FamilyCount		= "FAMILY_COUNT";				///<Number of families. This is actually the same as individual count
char *Config::ConfigKeywords::DatasetCount		= "DATASET_COUNT";
char *Config::ConfigKeywords::AddDataset		= "DATASET";
char *Config::ConfigKeywords::AddFamily			= "PED";
char *Config::ConfigKeywords::AddFamilyRandom	= "FLEXIPED";
char *Config::ConfigKeywords::AddCC				= "CC";

/**
 * This is used to change the pedigree header from 10 to 6 columns
 */
char *Config::ConfigKeywords::StandardPedigreeHeader = "USE_STD_PEDIGREE_HEADER";

/////Chromosom/block definition
char *Config::ConfigKeywords::SetDefaultBlock 	= "DEFAULT_BLOCK";			///<Sets up the default block
char *Config::ConfigKeywords::AddChromosome 	= "ADD_CHROMOSOME";				///<Start the definition of a new chromosome
char *Config::ConfigKeywords::LoadChromosome	= "LOAD_CHROMOSOME";			///<Chromosome loci are defined in file
char *Config::ConfigKeywords::AddBlock			= "ADD_BLOCK";					///<Add a block to the current chromosome

char *Config::ConfigKeywords::FirstDropPoint 	= "FIRST_DROP_POINT";			///<Specify where you want to start dumping the gene pools
char *Config::ConfigKeywords::DropFrequency 	= "DROP_FREQUENCY";				///<Specify how many generations are processed between the drops
char *Config::ConfigKeywords::DropCount 		= "DROP_COUNT";					///<How many drops are made after the initial drop is made
char *Config::ConfigKeywords::RandomizeAlleleFreq = "RANDOMIZE_ALLELE_FREQ";	///<Turn on/off allele frequency randomization


Config::Config()  {
	growthModel = new GrowthRate;
}


Config::~Config()
{
	delete growthModel;
}

PoolManager *Config::GetPoolManager() {
	return &pools;
}

GrowthRate *Config::GetGrowthRate() {
	return growthModel;
}



bool Config::SetValues(const char *key, const char *value, const char *line) {
	bool success=true;
	assert(strlen(key)>0);

	if 	    (strcmp(key, ConfigKeywords::PoolSize) == 0) 
		generalSettings.poolSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::DefaultAlleleFreq) == 0)
		SetDefaultAlleleFreq(line);
	else if (strcmp(key, ConfigKeywords::SetAlleleFrequency) == 0)
		SetAlleleFrequency(line);
	else if (strcmp(key, ConfigKeywords::SetAlleleFreqRange) == 0)
		SetAlleleFreqRange(line);
	else if (strcmp(key, ConfigKeywords::DatasetCount) == 0)
		datasetSettings.simSets = atoi(value);
	else if (strcmp(key, ConfigKeywords::RandomSeed) == 0) {
		generalSettings.seed = atoi(value);
		GAModel::seed = generalSettings.seed;
		Utility::Random::globalGenerator.Seed(generalSettings.seed);
	}
	else if (strcmp(key, ConfigKeywords::FamilyCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::FirstDropPoint) == 0)
		generalSettings.firstDropPoint = atoi(value);
	else if (strcmp(key, ConfigKeywords::DropFrequency) == 0)
		generalSettings.dropFrequency = atoi(value);
	else if (strcmp(key, ConfigKeywords::DropCount) == 0)
		generalSettings.dropCount = atoi(value);
//	else if (strcmp(key, ConfigKeywords::GenotypeError) == 0)
//		statusSettings.genotypeError = atof(value);
//	else if (strcmp(key, ConfigKeywords::PhenoCopy) == 0)
//		Simulation::DiseaseModel::phenocopyError = atof(value);
	else if (strcmp(key, ConfigKeywords::IndividualCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::AddModel) == 0)
		statusSettings.AddModel(line);
	else if (strcmp(key, ConfigKeywords::DefineModel) ==0)
		statusSettings.DefineModel(line, &pools);
	else if (strcmp(key, ConfigKeywords::GenerateModels) == 0)
		GAModel::doRunGA = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::AddChromosome) == 0)
		AddChromosome(line);
	else if (strcmp(key, ConfigKeywords::LoadChromosome) == 0)
		LoadChromosome(line);
	else if (strcmp(key, ConfigKeywords::AddBlock) == 0)
		AddBlock(line);
	else if (strcmp(key, ConfigKeywords::SetDefaultBlock) == 0)
		SetDefaultBlock(line);
	else if (strcmp(key, ConfigKeywords::RandomizeAlleleFreq) == 0)
		ChromPool::randomizeAlleleFreq = GetBoolean(value);
//	else if (strcmp(key, ConfigKeywords::AlleleLimits) == 0)
//		SetAlleleFrequency(line);
	else if (strcmp(key, ConfigKeywords::MaxPoolSize) == 0)
		GrowthRate::maxPoolSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::MinPoolSize) == 0)
		GrowthRate::minPoolSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::GrowthRate) == 0)
		SetGrowthRate(line);
	else if (strcmp(key, ConfigKeywords::AddDataset) == 0)
		AddDataset(line);
	else if (strcmp(key, ConfigKeywords::PathToJava) == 0)
		PoolManager::pathToJava=GetRemainderString(line);
	else if (strcmp(key, ConfigKeywords::PathToHaploview) ==0)
		PoolManager::pathToHaploview = GetRemainderString(line);
	else if (strcmp(key, ConfigKeywords::HaploviewSettings) == 0)
		PoolManager::haploviewSettings = GetRemainderString(line);
	else if (strcmp(key, ConfigKeywords::JavaSettings) == 0)
		PoolManager::javaSettings = GetRemainderString(line);
	else if (strcmp(key, ConfigKeywords::HaploviewOverviewSettings) == 0)
		PoolManager::haploviewOverviewSettings = GetRemainderString(line);
	else if (strcmp(key, ConfigKeywords::GenerateOverviewLD) == 0)
		PoolManager::generateOverviewLD = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::ProduceLDMaps) == 0)
		generalSettings.createLDMaps = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::KeepPhasedOutput) == 0)
		generalSettings.keepPhasedOutput = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::HaploviewWindowSize) == 0)
		PoolManager::haploWindowSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::HaploviewWindowStride) == 0)
		PoolManager::haploWindowStride = atoi(value);
	else if (strcmp(key, ConfigKeywords::GAFitnessThreshold) == 0) 
		GAModel::fitnessThreshold = atof(value);
	else if (strcmp(key, ConfigKeywords::GAIterations) == 0)
		GAModel::tries = atoi(value);
	else if (strcmp(key, ConfigKeywords::GenerateSampleLDMaps) ==0)
		generalSettings.createSampleLDmaps = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::MappingFN) == 0)
		SetMappingFn(value);
	else if (strcmp(key, ConfigKeywords::StandardPedigreeHeader) ==0)
		Individual::StandardPedigreeHeader = GetBoolean(value);
	else {
		cout<<"!! Unhandled line: "<<line<<"\n";
		success = false;
	}

	return success;
}

void Config::SetMappingFn(const char *val) {
	if (strcmp(val, "HALDANE") == 0){
		if (ChromPool::mappingFn)
			delete ChromPool::mappingFn;
		ChromPool::mappingFn = new HaldaneMapping();
	}
	else if (strcmp(val, "KOSAMBI") == 0) {
		if (ChromPool::mappingFn)
			delete ChromPool::mappingFn;
		ChromPool::mappingFn = new KosambiMapping();
	}
	else
		cout<<"Unknown mapping function: "<<val<<"\n";
}

string Config::GetRemainderString(const char *line) {
	string junk;
	string var;
	
	stringstream ss(line);
	ss>>junk;

	while (!ss.eof())  {
		ss>>junk;
		var += " " + junk;
		junk = "";
	}
	
	return var;

}
	

void Config::SetAlleleFreqRange(const char *line) {
	uint locus, chromID;
	float min, max;

	string cmd;
	stringstream ss(line);

	uint wordCount = CountColumns(line);
	if (wordCount == 5) {
		ss>>cmd>>chromID>>locus>>min>>max;
		if (!pools.ForceAlleleFreqRange(chromID - 1, locus - 1, min, max)) {
			cout<<"The requested chromosome ("<<chromID<<") hasn't been defined yet. Please define all chromosomes before setting allele frequencies\n";
			abort();
		}	
	}
	else  {
		cout<<"Unexpected format: "<<line<<"\n";
		cout<<"-- Unable to set allele frequency range\n";	

		cout<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tALLELE_FREQ_RANGE chromosome locus al1_min_freq al2_max_freq\n";
		abort();
	}
}	

void Config::SetAlleleFrequency(const char *line) {
	uint locus, chromID;
	float al1, al2;

	string cmd;
	stringstream ss(line);

	uint wordCount = CountColumns(line);
	if (wordCount == 5) {
		ss>>cmd>>chromID>>locus>>al1>>al2;
		if (!pools.ForceAlleleFrequency(chromID - 1, locus - 1, al1, al2)) {
			cout<<"The requested chromosome ("<<chromID<<") hasn't been defined yet. Please define all chromosomes before setting allele frequencies\n";
			abort();
		}	
	}
	else  {
		cout<<"Unexpected format: "<<line<<"\n";
		cout<<"-- Unable to set allele frequency\n";	

		cout<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tALLELE_FREQ chromosome locus freq1 freq2\n";
		abort();
	}
}

void Config::AddDataset(const char *line) {
	string cmd, type;
	stringstream ss(line);

	uint wordCount = CountColumns( line );

	ss>>cmd>>type;
	Sample *sample=NULL;

	if      (strcmp(type.c_str(), ConfigKeywords::AddFamily)==0) {
		uint aff=0, uaff=0;
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		ss>>aff>>uaff>>genotypeError>>phenoError>>missingData;
		
		if (aff == 0) {
			cout<<"No affected individuals defined. \n\t"<<line<<"\n";
			cout<<"Correct Syntax: \t\nDATASET PED affected_sibs unaffected_sibs genotype_error phenotype_error missing_data\n";
			abort();
		}
		
		sample=new PedigreeSample(aff, uaff, genotypeError, phenoError, missingData, 2);
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddFamilyRandom) == 0){ 
		uint maxFamSize=0;
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		ss>>maxFamSize>>genotypeError>>phenoError>>missingData;
		
		if (maxFamSize == 0) {
			cout<<"No upper limit to children count specified on FLEXIPED. \n\t"<<line<<"\n";
			cout<<"Correct Syntax: \t\nDATASET FLEXIPED affected_sibs unaffected_sibs genotype_error phenotype_error missing_data\n";
			abort();
		}
		
		sample=new PedigreeMixedSample(maxFamSize, genotypeError, phenoError, missingData);
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddCC) == 0) {
		float  percAffected=0.0, genotypeError=0.0, phenocopyError=0.0, missingData=0.0;
		ss>>percAffected>>genotypeError>>phenocopyError>>missingData;

		if (percAffected == 0.0) {
			cout<<"No affected individuals defined. \n\t"<<line<<"\n";
			abort();
		}

		sample=new BasicSample(percAffected, genotypeError, phenocopyError, missingData);
	}
	else 
		cout<<"Unknown dataset type: "<<type
			<<". Possible types include: "<<ConfigKeywords::AddFamily
			<<" and "<<ConfigKeywords::AddCC<<endl;
	if (sample)
		datasetSettings.samples.push_back(sample);		
}

void Config::SetDefaultAlleleFreq(const char *line) {
	string cmd;

	uint wordCount = CountColumns(line);

	if (wordCount != 3) {
		cout<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tDEFAULT_ALLELE_FREQ block_count min_recomb max_recomb\n";

		abort();
	}

	stringstream ss(line);
	ss>>cmd>>ChromPool::defFre1>>ChromPool::defFre2;
	
}

void Config::SetGrowthRate(const char *line) {
	string cmd;
	string type;

	stringstream ss(line);
	ss>>cmd;
	if (growthModel)
		delete growthModel;
	growthModel = GrowthRate::LoadModel(ss);

	cout<<"Writing growth table for "<<cmd<<" growth\n";
	string filename = generalSettings.outputName + string(".") + growthModel->GetType() + ".csv";
	ofstream file(filename.c_str(), ios_base::out);
	growthModel->DiagramGrowth(file, 1, 5000, 1);
}


void Config::LoadChromosome(const char *line) {
	uint wordCount = CountColumns(line);
	string filename, cmd;
	stringstream ss(line);
	ss>>cmd>>filename;
	currPool = pools.AddChromosome(filename.c_str());
	
}

void Config::AddChromosome(const char *line) {
	uint blockCount;
	float minR, maxR;
	string cmd;
	
	uint wordCount = CountColumns(line);

	if (wordCount != 4) {
		cout<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tADD_CHROMOSOME block_count min_recomb max_recomb\n";

		abort();
	}
	stringstream ss(line);
	ss>>cmd>>blockCount>>minR>>maxR;
	
	currPool = pools.AddChromosome(blockCount, minR, maxR, defaultBlock);
}


void Config::AddBlock(const char *line) {
	uint min, max;
	float blckMin, blckMax, snpMin, snpMax, prob;
	string cmd;

	uint wordCount = CountColumns(line);

	if (wordCount != 8) {
		cout<<"Misconfigured Block: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tADD_BLOCK min_loci_count max_loci_count block_recomb_min block_recomb_max snp_recomb_min snp_recomb_max probability\n";
		abort();
	}
	
	stringstream ss(line);
	ss>>cmd>>min>>max>>blckMin>>blckMax>>snpMin>>snpMax>>prob;

	currPool->DefineBlock(min, max, blckMin, blckMax, snpMin, snpMax, prob);
}

void Config::SetDefaultBlock(const char *line) {
	uint min, max;
	float blckMin=0.0, blckMax=0.0, snpMin=0.0, snpMax=0.0;
	stringstream ss(line);
	string cmd;

	uint wordCount = CountColumns(line);

	if (wordCount != 7) {
		cout<<"Misconfigured Block: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tDEFAULT_BLOCK min_loci_count max_loci_count block_recomb_min block_recomb_max snp_recomb_min snp_recomb_max\n";
		abort();
	}

	ss>>cmd>>min>>max>>blckMin>>blckMax>>snpMin>>snpMax;

	defaultBlock = ChromPool::BlockDefinition(min, max, blckMin, blckMax, snpMin, snpMax, 1.0, 0);
}


void Config::GenerateReport( ostream &os) {
	uint width = 45;
	
	os<<"\n\n-----------------------------------General Settings--------------------------\n";
	os<<setw(width)<<right<<"Random Seed: "<<generalSettings.seed<<endl;
	os<<setw(width)<<right<<"Pool Size: "<<generalSettings.poolSize<<endl;
	os<<setw(width)<<right<<"Starting Generation: "<<generalSettings.firstGeneration<<endl;
	if (generalSettings.firstGeneration > generalSettings.firstDropPoint)
		os<<setw(width)<<right<<"First Drop Point: "<<generalSettings.firstDropPoint<<endl;
	os<<setw(width)<<right<<"Drop Frequency: "<<generalSettings.dropFrequency<<endl;
	os<<setw(width)<<right<<"Drop Count: "<<generalSettings.dropCount<<endl;
	os<<setw(width + 8)<<right<<"Disease Models"<<endl;

	if (generalSettings.createLDMaps) {
		os<<setw(width)<<right<<"Production of LD maps using haploview: "<<"on\n";
		os<<setw(width)<<right<<"Path To Java: "<<PoolManager::pathToJava<<endl;
		os<<setw(width)<<right<<"Java settings: "<<PoolManager::javaSettings<<endl;
		os<<setw(width)<<right<<"Path To Haploview: "<<PoolManager::pathToHaploview<<endl;
		os<<setw(width)<<right<<"Haploview Settings: "<<PoolManager::haploviewSettings<<endl;
	}
	if (generalSettings.keepPhasedOutput)
		os<<setw(width)<<right<<"Phased chromosome files: "<<"Stored on filesystem\n";

	statusSettings.models.GenerateReport(os, width);

	os<<"\n\n-----------------------------------Pool Details------------------------------\n";
	pools.GenerateReport( os, width);
	
	os<<"\n\n-----------------------------Status Determination Settings-------------------\n";
	//float pe=Simulation::DiseaseModel::phenocopyError;
	//os<<setw(width)<<right<<"Phenocopy Error Rate: "<<pe<<endl;

	os<<"\n\n------------------------------Dataset Generation Settings--------------------\n";
	os<<setw(width)<<right<<"Datasets per Sample Type: "<<datasetSettings.simSets<<endl;
	os<<setw(width + 10)<<right<<"Dataset Samples"<<endl;
	for (uint i=0; i<datasetSettings.samples.size(); i++)
		datasetSettings.samples[i]->GenerateReport(os, width);
}


bool Config::PostLoad() {
	return pools.InitializePools(generalSettings.firstGeneration, growthModel->GetInitialPopulationSize(), generalSettings.outputName.c_str());
}

/**
 * @brief Check for missing data or impossible configuration settings
 */
bool Config::Validate() {
	bool isValid = true;

	

	return isValid;
}


}
