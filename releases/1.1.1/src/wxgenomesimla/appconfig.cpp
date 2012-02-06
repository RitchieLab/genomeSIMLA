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
#include "appconfig.h"
#include "utility/utility.h"
#include <assert.h>
#include <sstream>
#include "simulation/diseasemodel.h"
#include "simulation/pedigreesample.h"
#include "simulation/gamodel.h"
#include "simla/random.h"
#include "simulation/ldwriter.h"
#include "simulation/ldpngcomponent.h"
#include "simulation/ldplotter.h"
#include "utility/strings.h"
#include "simulation/lineargrowth.h"
#include <string>
#include <map>
namespace GenomeSIM {

using namespace std;

const char *AppConfig::AppConfigKeywords::GenerateOverviewLD 		= "GENERATE_LD_PLOTS";
const char *AppConfig::AppConfigKeywords::GAFitnessThreshold		= "GA_FITNESS_THRESHOLD";
const char *AppConfig::AppConfigKeywords::GAIterations				= "GA_ITERATIONS";
const char *AppConfig::AppConfigKeywords::MappingFN					= "SET_MAPPING_FN";
const char *AppConfig::AppConfigKeywords::MaxSnpsPerRow				= "MAX_SNPS_PER_ROW";
const char *AppConfig::AppConfigKeywords::MaxSNPDistance			= "MAX_SNP_DISTANCE";

////General Settings
///Number of unique expressions of each of the blocks to be found in the various pools
const char *AppConfig::AppConfigKeywords::PoolSize					= "POOL_SIZE";					
const char *AppConfig::AppConfigKeywords::MaxPoolSize				= "MAX_POOL_SIZE";				///<Ceiling for growth
const char *AppConfig::AppConfigKeywords::MinPoolSize				= "MIN_POOL_SIZE";				///<Floor
const char *AppConfig::AppConfigKeywords::RandomSeed				= "SEED";						///<Random seed
const char *AppConfig::AppConfigKeywords::AlleleLimits				= "ALLELE_LIMITS";				///<Range of allele limits
const char *AppConfig::AppConfigKeywords::GrowthRate				= "GROWTH_RATE";				///<Setup limits for growth rate / generation
const char *AppConfig::AppConfigKeywords::DefaultAlleleFreq			= "DEFAULT_ALLELE_FREQ";		///<The default frequencies of alleles
const char *AppConfig::AppConfigKeywords::SetAlleleFrequency 		= "ALLELE_FREQUENCY";			///<Force a frequency on a given locus
const char *AppConfig::AppConfigKeywords::SetAlleleFreqRange 		= "ALLELE_FREQ_RANGE";
const char *AppConfig::AppConfigKeywords::MinAlleleFreqThresh		= "MIN_ALLELE_FREQ_THRESH";
const char *AppConfig::AppConfigKeywords::BlockReportSize			= "BLOCK_REPORT_SIZE";			///<Number of blocks reported on

////Status Settings
const char *AppConfig::AppConfigKeywords::GenotypeError				= "GENOTYPE_ERROR";
const char *AppConfig::AppConfigKeywords::PhenoCopy 				= "PHENOCOPY_ERROR";
const char *AppConfig::AppConfigKeywords::SetModel					= "SET_MODEL";
const char *AppConfig::AppConfigKeywords::DefineModel				= "DEFINE_MODEL";
const char *AppConfig::AppConfigKeywords::LabelBased 				= "LABEL";
const char *AppConfig::AppConfigKeywords::IndexBased				= "INDEX";
const char *AppConfig::AppConfigKeywords::Simpen					= "SIMPEN";
const char *AppConfig::AppConfigKeywords::Simla						= "SIMLA";
const char *AppConfig::AppConfigKeywords::PenTable					= "PENTABLE";
const char *AppConfig::AppConfigKeywords::GenerateModels			= "GENERATE_MODELS";
const char *AppConfig::AppConfigKeywords::OffspringPerMating		= "OFFSPRING_PER_MATING";

/////Dataset settings
/**Number of individuals per file (for Case/control)*/
const char *AppConfig::AppConfigKeywords::IndividualCount			= "INDIVIDUAL_COUNT";	
/**Number of families. This is actually the same as individual count	*/
const char *AppConfig::AppConfigKeywords::FamilyCount				= "FAMILY_COUNT";				
const char *AppConfig::AppConfigKeywords::DatasetCount				= "DATASET_COUNT";
const char *AppConfig::AppConfigKeywords::AddDataset				= "DATASET";
const char *AppConfig::AppConfigKeywords::AddFamily					= "PED";
const char *AppConfig::AppConfigKeywords::AddFamilyRandom			= "FLEXIPED";
const char *AppConfig::AppConfigKeywords::AddFamilyType				= "FAMTYPE";
const char *AppConfig::AppConfigKeywords::AddCC						= "CC";
const char *AppConfig::AppConfigKeywords::IncludeAppConfiguration = "INCLUDE";
/**
 * This is used to change the pedigree header from 10 to 6 columns
 */
const char *AppConfig::AppConfigKeywords::StandardPedigreeHeader 	= "USE_STD_PEDIGREE_HEADER";

/////Chromosom/block definition
const char *AppConfig::AppConfigKeywords::SetDefaultBlock 	= "DEFAULT_BLOCK";			///<Sets up the default block
const char *AppConfig::AppConfigKeywords::AddChromosome 	= "ADD_CHROMOSOME";				///<Start the definition of a new chromosome
const char *AppConfig::AppConfigKeywords::LoadChromosome	= "LOAD_CHROMOSOME";			///<Chromosome loci are defined in file
const char *AppConfig::AppConfigKeywords::SeedChromosome	= "SEED_CHROMOSOME";
const char *AppConfig::AppConfigKeywords::AddBlock			= "ADD_BLOCK";					///<Add a block to the current chromosome

const char *AppConfig::AppConfigKeywords::FirstDropPoint 	= "FIRST_DROP_POINT";			///<Specify where you want to start dumping the gene pools
const char *AppConfig::AppConfigKeywords::DropFrequency 	= "DROP_FREQUENCY";				///<Specify how many generations are processed between the drops
const char *AppConfig::AppConfigKeywords::DropCount 		= "DROP_COUNT";					///<How many drops are made after the initial drop is made
const char *AppConfig::AppConfigKeywords::RandomizeAlleleFreq = "RANDOMIZE_ALLELE_FREQ";	///<Turn on/off allele frequency randomization
const char *AppConfig::AppConfigKeywords::DoDumpGenerationZero = "DUMP_GENERATION_ZERO";

const char *AppConfig::AppConfigKeywords::TrueTypeFontName		= "FONT";					///<The font to be used
const char *AppConfig::AppConfigKeywords::BlockReportBuffer		= "LD_REPORT_BUFFER_SIZE";	///<The number of snps to be used on either side of a block in an LD plot
const char *AppConfig::AppConfigKeywords::CssFilename			= "CSS_FILENAME";			///<Override default stylesheet
const char *AppConfig::AppConfigKeywords::ClosePoolsBetweenAdvancements = "CLOSE_POOLS_BETWEEN_DROPS";
const char *AppConfig::AppConfigKeywords::UseOriginalCrossing	= "USE_ORIGINAL_CROSSING";	///<This is for debugging purposes only (and benchmarking)
const char *AppConfig::AppConfigKeywords::PedUseOriginalCross  = "PED_USE_ORIGINAL_CROSSING";
const char *AppConfig::AppConfigKeywords::PhasedPedigrees		= "PHASED_PEDIGREES";
const char *AppConfig::AppConfigKeywords::UseAltLD				= "USE_ALT_LD";

const char *AppConfig::AppConfigKeywords::DrawRSquared			= "DRAW_RSQUARED_PLOTS";
const char *AppConfig::AppConfigKeywords::DrawDPrime			= "DRAW_DPRIME_PLOTS";
const char *AppConfig::AppConfigKeywords::WriteLDReport			= "WRITE_LD_REPORT";

const char *AppConfig::AppConfigKeywords::MaxLdIndividualCount	= "FAST_LD_POOL_SIZE";
const char *AppConfig::AppConfigKeywords::PlotScanSize			= "FAST_LD_PLOT_SIZE";



const char *AppConfig::AppConfigKeywords::BinaryDatasets		= "BINARY_DATASETS";
const char *AppConfig::AppConfigKeywords::LocSelection			= "LOCUS_SELECTOR";
const char *AppConfig::AppConfigKeywords::AddRegion 			= "ADD_REGION";

const char *AppConfig::AppConfigKeywords::MaxThreadCount		= "MAX_THREAD_COUNT";
const char *AppConfig::AppConfigKeywords::OverideRecThreadCount = "OVERRIDE_REC_THREAD_COUNT";

const char *AppConfig::AppConfigKeywords::SimultaneousChroms	= "SIMULTANEOUS_CHROM";
const char *AppConfig::AppConfigKeywords::ThreadsPerChrom		= "THREADS_PER_CHROM";

const char *AppConfig::AppConfigKeywords::NoFastLD				= "NO_FAST_LD";
const char *AppConfig::AppConfigKeywords::TargetPopSize			= "TARGET_POP_SIZE";

const char *AppConfig::AppConfigKeywords::LociPerChromReported	= "MAX_LOCI_PER_CHROM_REPORTED";

AppConfig *AppConfig::primaryInstance = NULL;


AppConfig::AppConfig() : currPool(NULL), doDumpGenerationZero(true), growthModelAppConfigured(false), configurationFilename("")  {
	growthModel = new LinearGrowth(100000, 0.0, 0.0);
	if (primaryInstance == NULL)
		primaryInstance = this;
}


AppConfig::~AppConfig()
{
	delete growthModel;
	if (primaryInstance == this)
		primaryInstance = NULL;
}

PoolManager *AppConfig::GetPoolManager() {
	return &pools;
}

GrowthRate *AppConfig::GetGrowthRate() {
	return growthModel;
}

void AppConfig::SetRandomSeed( uint seed ) {
		generalSettings.seed = seed;
		Simulation::StatusModel::GAModel::seed = generalSettings.seed;
		Utility::Random::globalGenerator.Seed(generalSettings.seed);
		Simla::Random::setseed(seed);
}

uint AppConfig::GetRandomSeed() {
	return generalSettings.seed;
}

bool AppConfig::WriteConfiguration(const char *config) {
	ofstream file(config);

	if (! file.good()) {
		cout<<"An error occured when trying to open the configuration file, "<<config<<"\n";
		return false;
	}

	file<<"# This file was created by GenomeSIMLA based on the selections from a graphical user interface (GUI).\n";
	file<<"# If using this configuration from the command line, make sure that the paths specified in this \n";
	file<<"# file are appropriate for the machine to be used to run the application properly. \n";

	int targetPopSize = ChromPool::targetPop;
	file<<AppConfigKeywords::PoolSize<<"\t"<<generalSettings.poolSize<<"\n";
	file<<AppConfigKeywords::DefaultAlleleFreq<<"\t"<<ChromPool::defFre1<<"\t"<<ChromPool::defFre2<<"\n";
	file<<AppConfigKeywords::DatasetCount<<"\t"<<datasetSettings.simSets<<"\n";
	file<<AppConfigKeywords::RandomSeed<<"\t"<<generalSettings.seed<<"\n";
	file<<AppConfigKeywords::FirstDropPoint<<"\t"<<generalSettings.firstDropPoint<<"\n";
	file<<AppConfigKeywords::DropFrequency<<"\t"<<generalSettings.dropFrequency<<"\n";
	file<<AppConfigKeywords::DropCount<<"\t"<<generalSettings.dropCount<<"\n";
	file<<AppConfigKeywords::TargetPopSize<<"\t"<<ChromPool::targetPop<<"\n";
	file<<AppConfigKeywords::SimultaneousChroms<<"\t"<<PoolManager::simultaneousChrom<<"\n";
	file<<AppConfigKeywords::ThreadsPerChrom<<"\t"<<PoolManager::threadsPerChrom<<"\n";
	file<<AppConfigKeywords::MaxPoolSize<<"\t"<<GrowthRate::maxPoolSize<<"\n";
	file<<AppConfigKeywords::MinPoolSize<<"\t"<<GrowthRate::minPoolSize<<"\n";
	file<<growthModel->GenerateCfgString()<<"\n";
	file<<AppConfigKeywords::MaxSnpsPerRow<<"\t"<<Simulation::Visualization::ImageParameters::maxSnpsPerRow<<"\n";
	file<<AppConfigKeywords::MaxSNPDistance<<"\t"<<Simulation::Visualization::LdPlotter::maxSnpDistance<<"\n";
	if (PoolManager::generateOverviewLD)
		file<<AppConfigKeywords::GenerateOverviewLD<<"\tYES\n";
	else
		file<<AppConfigKeywords::GenerateOverviewLD<<"\tNO\n";
//	file<<AppConfigKeywords::MappingFN<<"\t"<<ChromPool::mappingFn->GetType()<<"\n";
	file<<AppConfigKeywords::LociPerChromReported<<"\t"<<LocusSelection::maxReportedEntries<<"\n";
	file<<AppConfigKeywords::StandardPedigreeHeader<<"\t"<<Individual::StandardPedigreeHeader<<"\n";
	file<<AppConfigKeywords::BlockReportSize<<"\t"<<ChromPool::numberOfBlocksToReport<<"\n";
	file<<AppConfigKeywords::TrueTypeFontName<<"\t\""<<Simulation::Visualization::ImageParameters::font<<"\"\n";
	file<<AppConfigKeywords::BlockReportBuffer<<"\t"<<ChromPool::highBufferSize<<"\n";
	file<<AppConfigKeywords::CssFilename<<"\t"<<ChromPool::cssFilename<<"\n";
	file<<AppConfigKeywords::MaxLdIndividualCount<<"\t"<<ChromPool::maxLDIndividuals<<"\n";
	file<<AppConfigKeywords::PlotScanSize<<"\t"<<ChromPool::plotScanSize<<"\n";

	if (!ChromPool::fastLD) 
		file<<AppConfigKeywords::NoFastLD<<"\n";

	if (ChromPool::writeDPrimePlots) 
		file<<AppConfigKeywords::DrawDPrime<<"\tYES\n";
	else
		file<<AppConfigKeywords::DrawDPrime<<"\tNO\n";
	
	if (ChromPool::writeRSquaredPlots) 
		file<<AppConfigKeywords::DrawRSquared<<"\tYES\n";
	else
		file<<AppConfigKeywords::DrawRSquared<<"\tNO\n";
	
	if (ChromPool::writeLdTextReport)
		file<<AppConfigKeywords::WriteLDReport<<"\tYES\n";
	else
		file<<AppConfigKeywords::WriteLDReport<<"\tNO\n";
	

	//Locus Creation
	file<<"\n///////////////// Locus Details ///////////////////////////////////\n";
	file<<AppConfigKeywords::SetDefaultBlock<<"\t"<<
		defaultBlock.minSnpCount<<" "<<
		defaultBlock.maxSnpCount<<" "<<
		defaultBlock.minBlckMap<<" "<<
		defaultBlock.maxBlckMap<<" "<<
		defaultBlock.minSnpMap<<" "<<
		defaultBlock.maxSnpMap<<"\n";

	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();

	if (pool == NULL)
		cout<<"Well, we have no chromosomes to write....\n";

	while (pool) {
		cout<<"Writing a chromosome\n";
		pool->WriteConfiguration(file);		
		pool = itr.GetNext();
	}



	//Locus Selection
	generalSettings.WriteLocusSelections(file);

	//Datasets
	datasetSettings.WriteDatasetSelections(file);

	//Disease Model
	statusSettings.WriteModelSelection(file);


	return true;
}

bool AppConfig::SetValues(const char *key, const char *value, const char *line) {
	bool success=true;
	assert(strlen(key)>0);

	if 	    (strcmp(key, AppConfigKeywords::PoolSize) == 0) 
		generalSettings.poolSize = atoi(value);
	else if (strcmp(key, AppConfigKeywords::DefaultAlleleFreq) == 0)
		SetDefaultAlleleFreq(line);
	else if (strcmp(key, AppConfigKeywords::SetAlleleFrequency) == 0)
		SetAlleleFrequency(line);
	else if (strcmp(key, AppConfigKeywords::SetAlleleFreqRange) == 0)
		SetAlleleFreqRange(line);
	else if (strcmp(key, AppConfigKeywords::DatasetCount) == 0)
		datasetSettings.simSets = atoi(value);
	else if (strcmp(key, AppConfigKeywords::RandomSeed) == 0) 
		SetRandomSeed(atoi(value));
	else if (strcmp(key, AppConfigKeywords::FamilyCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, AppConfigKeywords::FirstDropPoint) == 0)
		generalSettings.firstDropPoint = atoi(value);
	else if (strcmp(key, AppConfigKeywords::DropFrequency) == 0)
		generalSettings.dropFrequency = atoi(value);
	else if (strcmp(key, AppConfigKeywords::DropCount) == 0)
		generalSettings.dropCount = atoi(value);
	else if (strcmp(key, AppConfigKeywords::TargetPopSize) ==0)
		ChromPool::targetPop = atoi(value);
	else if (strcmp(key, AppConfigKeywords::IndividualCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, AppConfigKeywords::DefineModel) ==0)
		statusSettings.DefineModel(line);
	else if (strcmp(key, AppConfigKeywords::GenerateModels) == 0)
		GAModel::doRunGA = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::AddChromosome) == 0)
		AddChromosome(line);
	else if (strcmp(key, AppConfigKeywords::LoadChromosome) == 0)
		LoadChromosome(line);
	else if (strcmp(key, AppConfigKeywords::SeedChromosome) == 0)
		SeedChromosome(line);
	else if (strcmp(key, AppConfigKeywords::AddBlock) == 0)
		AddBlock(line);
	else if (strcmp(key, AppConfigKeywords::SetDefaultBlock) == 0)
		SetDefaultBlock(line);
	else if (strcmp(key, AppConfigKeywords::RandomizeAlleleFreq) == 0)
		ChromPool::randomizeAlleleFreq = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::SimultaneousChroms) == 0)
		PoolManager::simultaneousChrom = atoi(value);
	else if (strcmp(key, AppConfigKeywords::ThreadsPerChrom) == 0)
		PoolManager::threadsPerChrom = atoi(value);
	else if (strcmp(key, AppConfigKeywords::NoFastLD) == 0)
		ChromPool::fastLD = false;
	else if (strcmp(key, AppConfigKeywords::OffspringPerMating) == 0) {
		stringstream ss(line);
		string cmd;
		uint min, max;

		ss>>cmd>>min>>max;
		ChromPool::SetOffspringLimits( min, max );	
	}

	else if (strcmp(key, AppConfigKeywords::MaxPoolSize) == 0) {
		if (growthModelAppConfigured) {
			cerr<<"Attempting to set MAX_POOL_SIZE after the a model was defined (GROWTH_RATE). \nPlease make sure that the bounds are set prior to the growth rate.\n";
			abort();
		}
		GrowthRate::maxPoolSize = atoi(value);
	}
	else if (strcmp(key, AppConfigKeywords::MinPoolSize) == 0) {
		if (growthModelAppConfigured) {
			cerr<<"Attempting to set MIN_POOL_SIZE after the a model was defined (GROWTH_RATE). \nPlease make sure that the bounds are set prior to the growth rate.\n";
			abort();
		}

		GrowthRate::minPoolSize = atoi(value);
	}
	else if (strcmp(key, AppConfigKeywords::GrowthRate) == 0)
		SetGrowthRate(line);
	else if (strcmp(key, AppConfigKeywords::AddDataset) == 0)
		AddDataset(line);
	else if (strcmp(key, AppConfigKeywords::MaxSnpsPerRow) == 0)
		Simulation::Visualization::ImageParameters::maxSnpsPerRow = atoi(value);
	else if (strcmp(key, AppConfigKeywords::MaxSNPDistance) == 0)
		Simulation::Visualization::LdPlotter::maxSnpDistance = atoi(value);
	else if (strcmp(key, AppConfigKeywords::GenerateOverviewLD) == 0)
		PoolManager::generateOverviewLD = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::GAFitnessThreshold) == 0) 
		Simulation::StatusModel::GAModel::fitnessThreshold = atof(value);
	else if (strcmp(key, AppConfigKeywords::GAIterations) == 0)
		Simulation::StatusModel::GAModel::tries = atoi(value);
	else if (strcmp(key, AppConfigKeywords::MappingFN) == 0)
		SetMappingFn(value);
	else if (strcmp(key, AppConfigKeywords::LociPerChromReported) == 0)
		LocusSelection::maxReportedEntries = atoi(value);
	else if (strcmp(key, AppConfigKeywords::StandardPedigreeHeader) ==0)
		Individual::StandardPedigreeHeader = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::DoDumpGenerationZero) == 0)
		doDumpGenerationZero = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::BlockReportSize) == 0) 
		ChromPool::numberOfBlocksToReport = atoi(value);
	else if (strcmp(key, AppConfigKeywords::TrueTypeFontName) ==0) {
		stringstream s(line);
		string junk;
		s>>junk;
		Simulation::Visualization::ImageParameters::font = ParseFilename(s, "Font Filename");
	}
	else if (strcmp(key, AppConfigKeywords::BlockReportBuffer) ==0)
		ChromPool::highBufferSize = atoi(value);
	else if (strcmp(key, AppConfigKeywords::CssFilename) == 0)
		ChromPool::cssFilename = value;
	else if (strcmp(key, AppConfigKeywords::UseOriginalCrossing) == 0)
		ChromPool::UseOriginalCrossing = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::PedUseOriginalCross) == 0)
		PedigreeSample::UseOriginalCross = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::PhasedPedigrees) == 0)
		Individual::PhasedPedigrees = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::UseAltLD) == 0)
		ChromPool::UseAltLD = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::MaxLdIndividualCount) == 0)
		ChromPool::maxLDIndividuals = atoi(value);
	else if (strcmp(key, AppConfigKeywords::PlotScanSize) == 0)
		ChromPool::plotScanSize = atoi(value);
	else if (strcmp(key, AppConfigKeywords::DrawDPrime) == 0)
		ChromPool::writeDPrimePlots = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::DrawRSquared) == 0)
		ChromPool::writeRSquaredPlots = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::WriteLDReport) == 0)
		ChromPool::writeLdTextReport  = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::ClosePoolsBetweenAdvancements) == 0)
		generalSettings.closePoolsBetweenAdvancing = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::BinaryDatasets) == 0)
		datasetSettings.binaryDatasets = GetBoolean(value);
	else if (strcmp(key, AppConfigKeywords::MinAlleleFreqThresh) == 0) {
		Individual::minAlFreqThreshold = atoi(value);
		if (Individual::minAlFreqThreshold > 0.0) {
			cout<<"!!!!!!!!!!!!!!!!!!!!!\n!! MIN_ALLELE_FREQ_THRESH is currently not recommended when generating datasets. Be warned\n";
		}
	}
	else if (strcmp(key, AppConfigKeywords::IncludeAppConfiguration) == 0) {
		if (!LoadSettings(value))  {
			cerr<<"Unable to load settings from the file: "<<value<<". Refusing to continue\n";
			abort();
		}	
	}
	else if (strcmp(key, AppConfigKeywords::LocSelection) == 0) 
		generalSettings.AddSelector( line );
	else if (strcmp(key, AppConfigKeywords::AddRegion) == 0) 
		generalSettings.AddSelectorRegion( line );
	else {
		cerr<<"!! Unhandled line: "<<line<<"\n";
		success = false;
	}

	return success;
}

void AppConfig::SetMappingFn(const char *val) {
/*	if (strcmp(val, "HALDANE") == 0){
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
*/
}

string AppConfig::GetRemainderString(const char *line) {
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
	

void AppConfig::SetAlleleFreqRange(const char *line) {
	uint locus, chromID;
	float min, max;

	string cmd;
	stringstream ss(line);

	uint wordCount = CountColumns(line);
	if (wordCount == 5) {
		ss>>cmd>>chromID>>locus>>min>>max;
		if (!pools.ForceAlleleFreqRange(chromID - 1, locus - 1, min, max)) {
			cerr<<"The requested chromosome ("<<chromID<<") hasn't been defined yet. Please define all chromosomes before setting allele frequencies\n";
			abort();
		}	
	}
	else  {
		cerr<<"Unexpected format: "<<line<<"\n";
		cerr<<"-- Unable to set allele frequency range\n";	

		cerr<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cerr<<"Correct syntax: \n\tALLELE_FREQ_RANGE chromosome locus al1_min_freq al2_max_freq\n";
		abort();
	}
}	

void AppConfig::SetAlleleFrequency(const char *line) {
	uint locus, chromID;
	float al1, al2;

	string cmd;
	stringstream ss(line);

	uint wordCount = CountColumns(line);
	if (wordCount == 5) {
		ss>>cmd>>chromID>>locus>>al1>>al2;
		if (!pools.ForceAlleleFrequency(chromID - 1, locus - 1, al1, al2)) {
			cerr<<"The requested chromosome ("<<chromID<<") hasn't been defined yet. Please define all chromosomes before setting allele frequencies\n";
			abort();
		}	
	}
	else  {
		cerr<<"Unexpected format: "<<line<<"\n";
		cerr<<"-- Unable to set allele frequency\n";	

		cerr<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cerr<<"Correct syntax: \n\tALLELE_FREQ chromosome locus freq(A) freq(a)\n";
		abort();
	}
}
void AppConfig::ClearDatasets() {
	datasetSettings.Purge();
}
void AppConfig::AddDataset(const char *line) {
	string cmd, type;
	uint columns = CountColumns(line);
	stringstream ss(line);

	//uint wordCount = CountColumns( line );

	ss>>cmd>>type;
	Sample *sample=NULL;

	static PedigreeSample* lastPedigree = NULL;			///<Help keep up with the last sample added

	if (strcmp(type.c_str(), AppConfigKeywords::AddFamily)==0) {
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		string desc;
		ss>>desc>>genotypeError>>phenoError>>missingData;
		
		
		lastPedigree=new PedigreeSample(genotypeError, phenoError, missingData, desc.c_str());
		sample = lastPedigree;
	}
	else if (strcmp(type.c_str(), AppConfigKeywords::AddFamilyType) == 0) {
		uint aff=0, uaff=0, rnd=0, count=0;
		ss>>aff>>uaff>>rnd>>count;
		if (lastPedigree) {
			lastPedigree->AddFamilyType(aff, uaff, count,  rnd);

			if (aff + uaff == 0 || count == 0) {
				cerr<<"Invalid Family Type Definition:\n\t"<<line<<"\n";
				cerr<<"Correct Syntax: \t\nDATASET PED affected_sibs unaffected_sibs extra_sibs genotype_error phenotype_error missing_data\n";
				abort();
			}
		}
		else {
			cerr<<"Error: "<<cmd<<" "<<type<<" was used before a family was defined\n";
			abort();
		}

	}
	else if (strcmp(type.c_str(), AppConfigKeywords::AddFamilyRandom) == 0){ 
		uint maxFamSize=0;
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		string desc;

		ss>>desc>>maxFamSize>>genotypeError>>phenoError>>missingData;
		
		if (maxFamSize == 0) {
			cerr<<"No upper limit to children count specified on FLEXIPED. \n\t"<<line<<"\n";
			cerr<<"Correct Syntax: \t\nDATASET FLEXIPED affected_sibs unaffected_sibs genotype_error phenotype_error missing_data\n";
			abort();
		}
		
		sample=new PedigreeMixedSample(maxFamSize, genotypeError, phenoError, missingData, desc.c_str());
	}
	else if (strcmp(type.c_str(), AppConfigKeywords::AddCC) == 0) {
		uint affCount = 0, unaffCount = 0;
		float  genotypeError=0.0, phenocopyError=0.0, missingData=0.0;
		string desc;

		if (columns != 8) {
			cerr<<"Unexpected format: "<<line<<"\n";
			cerr<<"-- Unable to configure case/control\n";	
	
			cerr<<"Correct syntax: \n\tDATASET CC name [num aff] [num unaff] [genotype error] [phenocopy error] [missing data]\n";
			abort();
		}
		ss>>desc>>affCount>>unaffCount;
		ss>>genotypeError>>phenocopyError>>missingData;
	
		if (affCount < 1)
			cerr<<"!! Due to the absence of affecteds present, no status column will be written\n";

		sample=new BasicSample(affCount, unaffCount, 0, genotypeError, phenocopyError, missingData, desc.c_str());
	}
	else 
		cout<<"Unknown dataset type: "<<type
			<<". Possible types include: "<<AppConfigKeywords::AddFamily
			<<" and "<<AppConfigKeywords::AddCC<<endl;
	if (sample)
		datasetSettings.samples.push_back(sample);		
}

void AppConfig::SetDefaultAlleleFreq(const char *line) {
	string cmd;

	uint wordCount = CountColumns(line);

	if (wordCount != 3) {
		cerr<<"Misconfigured default allele frequency: \n\t"<<line<<"\n";
		cerr<<"Correct syntax: \n\tDEFAULT_ALLELE_FREQ min_freq max_freq\n";

		abort();
	}

	stringstream ss(line);
	ss>>cmd>>ChromPool::defFre1>>ChromPool::defFre2;
	
}

void AppConfig::SetGrowthRate(const char *line) {
	string cmd;
	string type;

	stringstream ss(line);
	ss>>cmd;
	if (growthModel)
		delete growthModel;
	growthModel = GrowthRate::LoadModel(ss);
	growthModelAppConfigured = true;
/*	cout<<"Writing growth table for "<<cmd<<" growth\n";
	string filename = generalSettings.outputName + string(".") + growthModel->GetType() + ".csv";
	ofstream file(filename.c_str(), ios_base::out);
	growthModel->DiagramGrowth(file, 1, 5000, 1);
	*/
}


bool AppConfig::LoadSettings(const char *filename) {
	configurationFilename = filename;
	char c  = '#';
	Utility::LineParser lp(c);
	return lp.Parse(filename, this) > 0;
}


//SEED_CHROMOSOME chr22.ped chr22.loc Chrom-22
void AppConfig::SeedChromosome(const char *line) {
	string cmd, poolDataFilename, fmt, locFilename, label;
	int headerCols = 0;

	uint wordCount = CountColumns(line);

	if (wordCount < 5 || wordCount > 6) {
		cerr<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		cerr<<"Correct Syntax: \n\tSEED_CHROMOSOME PED/CC NumHeaderCols source_file loc_file [label]\n";
		cerr<<"\tExample: SEED_CHROMOSOME PED 6 chr22.ped chr22.loc Chrom-22\n";
		abort();
	}

	stringstream ss(line);

	ss>>cmd>>fmt>>headerCols>>poolDataFilename>>locFilename;
	
	if (wordCount > 3)
		ss>>label;
	else
		label = Utility::ExtractBaseFilename( poolDataFilename.c_str() );

	currPool = pools.SeedChromosome(poolDataFilename.c_str(), locFilename.c_str(), label.c_str(), fmt.c_str(), headerCols);

}

void AppConfig::LoadChromosome(const char *line) {
	string filename, cmd, label;
	stringstream ss(line);
	
	ss>>cmd;
	filename = ParseFilename(ss, "Chromosome File");

	ss>>label;
	while (!ss.eof()) {
		string s;
		ss>>s;
		label=label + " " + s;
	}

	if (filename.length() < 1) {
		stringstream ss;
		ss<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		ss<<"Correct syntax: \n\tLOAD_CHROMOSOME chromosome_file [label]\n";
		
		throw Utility::Exception::General(ss.str().c_str());
	}

	if (label.length() < 1)
		label = Utility::ExtractBaseFilename( filename.c_str());
	
	currPool = pools.AddChromosome(filename.c_str(), label.c_str());

}

void AppConfig::AddChromosome(const char *line) {
	uint blockCount;
	string cmd;
	string label;
	
	uint wordCount = CountColumns(line);

	if (wordCount < 2 || wordCount > 4) {
		cerr<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		cerr<<"Correct syntax: \n\tADD_CHROMOSOME block_count [label]\n";

		abort();
	}
	stringstream ss(line);
	ss>>cmd>>blockCount;

	if (wordCount == 3) {
		ss>>label;
	}
	else {
		char *lbl = new char[64];
		sprintf(lbl, "Chromosome-%d", (int)pools.GetPoolCount());
		label = lbl;
		delete[] lbl;
	}
	
	currPool = pools.AddChromosome(blockCount, defaultBlock, label.c_str());
	
}	


void AppConfig::AddBlock(const char *line) {
	if (currPool) {
		uint min, max;
		float blckMin, blckMax, snpMin, snpMax, prob;
		string cmd;
	
		uint wordCount = CountColumns(line);
	
		if (wordCount != 8) {
			cerr<<"Misconfigured Block: \n\t"<<line<<"\n";
			cerr<<"Correct syntax: \n\tADD_BLOCK min_loci_count max_loci_count block_recomb_min block_recomb_max snp_recomb_min snp_recomb_max probability\n";
			abort();
		}
		
		stringstream ss(line);
		ss>>cmd>>min>>max>>blckMin>>blckMax>>snpMin>>snpMax>>prob;
	
		currPool->DefineBlock(min, max, blckMin, blckMax, snpMin, snpMax, prob);
	}
	else {
		cerr<<"\nError. A Chromosome must be defined before block definitions are added to it\n";
		abort();
	}
}

void AppConfig::SetDefaultBlock(ChromPool::BlockDefinition& defBlock) {
	defaultBlock = defBlock;
}
void AppConfig::SetDefaultBlock(const char *line) {
	uint min, max;
	float blckMin=0.0, blckMax=0.0, snpMin=0.0, snpMax=0.0;
	stringstream ss(line);
	string cmd;

	uint wordCount = CountColumns(line);

	if (wordCount != 7) {
		cerr<<"Misconfigured Block: \n\t"<<line<<"\n";
		cerr<<"Correct syntax: \n\tDEFAULT_BLOCK min_loci_count max_loci_count block_recomb_min block_recomb_max snp_recomb_min snp_recomb_max\n";
		abort();
	}

	ss>>cmd>>min>>max>>blckMin>>blckMax>>snpMin>>snpMax;

	defaultBlock = ChromPool::BlockDefinition(min, max, blckMin, blckMax, snpMin, snpMax, 1.0, 0);
}


void AppConfig::GenerateReport( ostream &os) {
	uint width = 45;
	
	os<<"\n\n-----------------------------------General Settings--------------------------\n";
	os<<setw(width)<<right<<"Random Seed: "<<generalSettings.seed<<endl;
	os<<setw(width)<<right<<"Starting Generation: "<<generalSettings.firstGeneration<<endl;
	if (generalSettings.firstGeneration < generalSettings.firstDropPoint)
		os<<setw(width)<<right<<"First Drop Point: "<<generalSettings.firstDropPoint<<endl;

	os<<setw(width)<<right<<"Drop Frequency: "<<generalSettings.dropFrequency<<endl;
	os<<setw(width)<<right<<"Drop Count: "<<generalSettings.dropCount<<endl;

	os<<setw(width)<<right<<"Growth Report: "<<GetGrowthChartFilename()<<"\n";


	os<<"\n";
	os<<"\n\n-----------------------------------Graphical Reports--------------------------\n";
	os<<setw(width)<<right<<"Production of LD maps: ";
	if (PoolManager::generateOverviewLD) 	{
		os<<"On\n";
		ChromPool::GenerateReport(os, width + 5);
	}
	else
		os<<"Off\n";




	os<<"\n\n------------------------------------Threading--------------------------------\n";
	os<<setw(width)<<right<<"Simultaneous Chromosomes: "<<PoolManager::simultaneousChrom<<"\n";
	os<<setw(width)<<right<<"Threads Per Chromosome: "<<PoolManager::threadsPerChrom<<"\n\n";

	os<<"\n\n-----------------------------------Pool Details------------------------------\n";
	pools.GenerateReport( os, width);
	
	os<<"\n\n-----------------------------Status Determination Settings-------------------\n";
	if (statusSettings.model) 
		//statusSettings.models.GenerateReport(os, width);
		statusSettings.model->GenerateReport(os, width);
	else
		cout<<setw(width)<<right<<"No model selected\n";
	//float pe=Simulation::DiseaseModel::phenocopyError;
	//os<<setw(width)<<right<<"Phenocopy Error Rate: "<<pe<<endl;

	os<<"\n\n------------------------------Dataset Generation Settings--------------------\n";
	os<<setw(width)<<right<<"Datasets per Sample Type: "<<datasetSettings.simSets<<endl;
	os<<setw(width + 10)<<right<<"Dataset Samples"<<endl;
	for (uint i=0; i<datasetSettings.samples.size(); i++)
		datasetSettings.samples[i]->GenerateReport(os, width);
}


void AppConfig::SummarizeDiseaseModel(ostream &os, vector<Locus *>diseaseLoci) {
/*	uint estimatedGenerations = generalSettings.firstDropPoint + generalSettings.dropCount * generalSettings.dropFrequency;
	bool success = pools.InitializeLoci(generalSettings.firstGeneration, generalSettings.outputName.c_str(), generalSettings.doLoad);
	if (success) {
		cout<<"."; cout.flush();
*/

	if (statusSettings.model) {
		statusSettings.model->Refresh(&pools);
		statusSettings.model->Load();
		statusSettings.model->GenerateDetailedReport(os, diseaseLoci);
	}
	else
		os<<"<B>No Model Defined</B>\n";
}

bool AppConfig::PostLoad() {
	cout<<"Intializing (this may take a few minutes)";
	uint estimatedGenerations = generalSettings.firstDropPoint + generalSettings.dropCount * generalSettings.dropFrequency;
	generalSettings.growthChartFilename = growthModel->DrawGrowthChart(generalSettings.outputName.c_str(), 1, estimatedGenerations, (uint)(estimatedGenerations *.01));
	bool success = pools.InitializeLoci(generalSettings.firstGeneration, generalSettings.outputName.c_str(), true);
	//bool success = pools.InitializeLoci(generalSettings.firstGeneration, generalSettings.outputName.c_str(), generalSettings.doLoad);
	if (success) {
		cout<<"."; cout.flush();
		PenetranceModel *model = statusSettings.InitModel(&pools, datasetSettings.writeDatasets);
		if (model) 
			pools.InsertModel(model);

		//Now we want to finalize the locus selection stuff 
		generalSettings.ReconcileSelections(pools);

		success =	pools.InitializePools(generalSettings.firstGeneration, growthModel->GetInitialPopulationSize(), generalSettings.outputName.c_str(), generalSettings.doLoad, doDumpGenerationZero || generalSettings.closePoolsBetweenAdvancing, generalSettings.closePoolsBetweenAdvancing);
		cout<<"."; cout.flush();
	}
	cout<<".Done!\n";
	return success;
}

/**
 * @brief Check for missing data or impossible configuration settings
 */
bool AppConfig::Validate() {
	bool isValid = true;

	

	return isValid;
}


PenetranceModel *AppConfig::Status::LoadModel(const char *project, PoolManager *poolMgr) {
	poolMgr->PrepForSampling(project, poolMgr->GetCurrentGeneration(), model);

	model->Refresh(poolMgr);
	model->Load();
	return model;
}


void AppConfig::Status::WriteModelSelection(ostream &file) {	
	if (model) {
		file<<"\n///////////////// Disease model ///////////////////////////////////\n";
			
		file<<model->GetModelConfiguration()<<"\n";
	}
}

PenetranceModel *AppConfig::Status::InitModel(PoolManager *poolMgr, bool requireModel) {
	stringstream ss(modelCfg.c_str());

	string cmd;
	
	ss>>cmd;			//>>percentage;

	if (model)  {
		cout<<"Overwriting previous model definition. (";
		model->GenerateReport(cout, 0);
		cout<<")\n";
		delete model;
	}

	//Right now, all models are based on penetranceModels- which is a factory
	model = Simulation::StatusModel::PenetranceModel::GetModel(ss, poolMgr);

	if (requireModel && model==NULL) {
		if (modelCfg=="")
			cerr<<"There was no model specified. Please set up a disease model in the configuration file\n";
		else
			cerr<<"Invalid model definition: "<<modelCfg<<". Does this file really exist?\n";
		abort();
	}
	return model;
	//model->Init(ss, poolMgr);
}

void AppConfig::ConfigureDiseaseModel(bool writeDatasets) {
	PenetranceModel *model = statusSettings.InitModel(&pools, writeDatasets);
	if (model) 
		pools.InsertModel(model);

	//Now we want to finalize the locus selection stuff 
	generalSettings.ReconcileSelections(pools);

}

void AppConfig::General::PartiallyReconcileSelections(PoolManager &pools) {
	LocusMap &map = pools.GetLocusMap();
	//work through each of the region lines and add them to the appropriate selector
	vector<string>::iterator itr = selectionRegions.begin();
	vector<string>::iterator end = selectionRegions.end();

	while (itr != end) {
		stringstream line(itr->c_str());
		cout<<"Trying to parse the following line: "<<itr->c_str()<<"\n";
		string cmd, regionLbl, label1, label2;
		line>>cmd>>regionLbl>>label1>>label2;

		if (regionLbl.length() > 0 && label1.length() > 0 && label2.length() > 0) {
			if (map.find(label1.c_str()) == map.end() || map.find(label2.c_str()) == map.end()) {
				cerr<<"No locus was found to match, "<<label1<<". Please correct this and run again. \n\t"<<line.str()<<"\n";
				abort();
			}
			if (locusSelections.find(regionLbl.c_str()) == locusSelections.end() ) {
				cerr<<"No matching selector was found, "<<regionLbl<<". Please correct this and run again. \n\t"<<line.str()<<"\n";
				abort();
			}
			Locus l1 = map[label1.c_str()];
			Locus l2 = map[label2.c_str()];
			locusSelections[regionLbl.c_str()].AddRange(l1, l2);
		}
		else
			cout<<"Malformed Region specification: "<<line.str()<<"\n";
		itr++;
	}

}

void AppConfig::General::ReconcileSelections(PoolManager &pools) {
	LocusMap &map = pools.GetLocusMap();
	//work through each of the region lines and add them to the appropriate selector
	vector<string>::iterator itr = selectionRegions.begin();
	vector<string>::iterator end = selectionRegions.end();

	while (itr != end) {
		stringstream line(itr->c_str());
		string cmd, regionLbl, label1, label2;
		line>>cmd>>regionLbl>>label1>>label2;

		if (regionLbl.length() > 0 && label1.length() > 0 && label2.length() > 0) {
			if (map.find(label1.c_str()) == map.end() || map.find(label2.c_str()) == map.end()) {
				cerr<<"No locus was found to match, "<<label1<<". Please correct this and run again. \n"<<line.str()<<"\n";
				abort();
			}
			if (locusSelections.find(regionLbl.c_str()) == locusSelections.end() ) {
				cerr<<"No matching selector was found, "<<regionLbl<<". Please correct this and run again. \n"<<line.str()<<"\n";
				abort();
			}
			Locus l1 = map[label1.c_str()];
			Locus l2 = map[label2.c_str()];
			locusSelections[regionLbl.c_str()].AddRange(l1, l2);
		}
		itr++;
	}

	Simulation::Visualization::LocusSelection defaultRegions = locusSelections["ALL_SELECTORS"];
	Simulation::LocusSelectionArray selArray;

	LocusSelectionMap::iterator cur = locusSelections.begin();
	LocusSelectionMap::iterator selEnd = locusSelections.end();

	while (cur != selEnd) {
		if (cur->first != "ALL_SELECTORS") {
			cur->second.ExtractRanges(defaultRegions);
			selArray.push_back(cur->second);
		}
		cur++;
	}

	//Finally, we can set the pool manager's locus selection objects appropriately
	pools.SetLocusSelectors(selArray);
}

void AppConfig::General::AddSelector(const char *selLine) {
	stringstream ss(selLine);

	string label;
	float mafTarget = 0.0, mafMin = 0.0, mafMax = 0.0;
	uint blTarget = 0, blMin = 0, blMax = 0;
	string command;	

	ss>>command>>label>>mafTarget>>mafMin>>mafMax>>blTarget>>blMin>>blMax;

	Simulation::Visualization::LocusSelection newSel(mafTarget, mafMin, mafMax, blTarget, blMin, blMax);

	if ((locusSelections.find(label) != locusSelections.end()) && label != "ALL_SELECTORS") {
		cerr<<"Locus Selector label, "<<label<<" must be unique.";
		abort();
	}	
	newSel.SetLabel(label.c_str());
	size_t lineLength = strlen(selLine);
	char *desc = new char[lineLength + 1];
	desc[lineLength] = '\n';
	ss.getline(desc, lineLength);
	newSel.SetDescription(desc);
	locusSelections[label] = newSel;

	delete[] desc;	
}

string AppConfig::Dataset::GetDatasetConfiguration() {
	static string yesno[] = { "NO", "YES" };
	stringstream file;
	vector<Sample*>::iterator itr = samples.begin();
	vector<Sample*>::iterator end = samples.end();

	while (itr != end) {
		Sample *sample = *itr;
		sample->WriteConfiguration(file);
		itr++;
	}
	file<<"DATASET_COUNT "<<simSets<<"\n";
	file<<"BINARY_DATASETS "<<yesno[binaryDatasets]<<"\n";
	file<<"USE_STD_PEDIGREE_HEADER "<<yesno[Individual::StandardPedigreeHeader]<<"\n";
	return file.str();
}

void AppConfig::Dataset::WriteDatasetSelections(ostream &file) {
	file<<"\n///////////////// Datasets //////////////////////////////////////\n";
	file<<GetDatasetConfiguration()<<"\n";
}

void AppConfig::General::WriteLocusSelections(ostream &os) {
	LocusSelectionMap::iterator itr = locusSelections.begin();
	LocusSelectionMap::iterator end = locusSelections.end();

	while (itr != end) {
		LocusSelection &sel = itr->second;

		if (itr->first != "ALL_SELECTORS") {
			os<<AppConfigKeywords::LocSelection<<" "<<sel.GetConfiguration(" ")<<"\n";
			size_t rangeCount = sel.GetRangeCount();
			for (size_t i=0; i<rangeCount; i++) {
				os<<AppConfigKeywords::AddRegion<<" "<<itr->first<<" ";
				SimpleRange<Locus> lRange = sel.GetRange(i);
				os<<lRange.min.GetLabel()<<" "<<lRange.max.GetLabel()<<"\n";
			}
		}
		itr++;
	}
}	


}
