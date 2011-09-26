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
#include "simulation/pedigreereferencesample.h"
#include "simulation/gamodel.h"
#include "simla/random.h"
#include "simulation/ldwriter.h"
#include "simulation/ldpngcomponent.h"
#include "simulation/ldplotter.h"
#include "utility/strings.h"
#include <string>
#include <map>
#include <sys/time.h>
#include "simulation/lineargrowth.h"
namespace GenomeSIM {

using namespace std;

const char *Config::ConfigKeywords::GenerateOverviewLD 	= "GENERATE_LD_PLOTS";
const char *Config::ConfigKeywords::GAFitnessThreshold	= "GA_FITNESS_THRESHOLD";
const char *Config::ConfigKeywords::GAIterations 			= "GA_ITERATIONS";
const char *Config::ConfigKeywords::MappingFN				= "SET_MAPPING_FN";
const char *Config::ConfigKeywords::MaxSnpsPerRow			= "MAX_SNPS_PER_ROW";
const char *Config::ConfigKeywords::MaxSNPDistance		= "MAX_SNP_DISTANCE";

////General Settings
///Number of unique expressions of each of the blocks to be found in the various pools
const char *Config::ConfigKeywords::PoolSize			= "POOL_SIZE";
const char *Config::ConfigKeywords::MaxPoolSize			= "MAX_POOL_SIZE";				///<Ceiling for growth
const char *Config::ConfigKeywords::MinPoolSize			= "MIN_POOL_SIZE";				///<Floor
const char *Config::ConfigKeywords::RandomSeed			= "SEED";						///<Random seed
const char *Config::ConfigKeywords::AlleleLimits		= "ALLELE_LIMITS";				///<Range of allele limits
const char *Config::ConfigKeywords::GrowthRate			= "GROWTH_RATE";				///<Setup limits for growth rate / generation
const char *Config::ConfigKeywords::DefaultAlleleFreq	= "DEFAULT_ALLELE_FREQ";		///<The default frequencies of alleles
const char *Config::ConfigKeywords::SetAlleleFrequency 	= "ALLELE_FREQUENCY";			///<Force a frequency on a given locus
const char *Config::ConfigKeywords::SetAlleleFreqRange 	= "ALLELE_FREQ_RANGE";
const char *Config::ConfigKeywords::MinAlleleFreqThresh = "MIN_ALLELE_FREQ_THRESH";
const char *Config::ConfigKeywords::BlockReportSize		= "BLOCK_REPORT_SIZE";			///<Number of blocks reported on

////Status Settings
const char *Config::ConfigKeywords::GenotypeError		= "GENOTYPE_ERROR";
const char *Config::ConfigKeywords::PhenoCopy 			= "PHENOCOPY_ERROR";
const char *Config::ConfigKeywords::SetModel			= "SET_MODEL";
const char *Config::ConfigKeywords::DefineModel			= "DEFINE_MODEL";
const char *Config::ConfigKeywords::LabelBased 			= "LABEL";
const char *Config::ConfigKeywords::IndexBased			= "INDEX";
const char *Config::ConfigKeywords::Simpen				= "SIMPEN";
const char *Config::ConfigKeywords::Simla				= "SIMLA";
const char *Config::ConfigKeywords::TemplateModel		= "TEMPLATE_MODEL";
const char *Config::ConfigKeywords::PenTable			= "PENTABLE";
const char *Config::ConfigKeywords::GenerateModels		= "GENERATE_MODELS";
const char *Config::ConfigKeywords::OffspringPerMating	= "OFFSPRING_PER_MATING";

const char *Config::ConfigKeywords::TemplateMatchAttempts	= "TEMPLATE_MATCH_COUNT";
const char *Config::ConfigKeywords::TemplateThresh		= "TEMPLATE_AFFECTED_THRESH";
/////Dataset settings
/**Number of individuals per file (for Case/control)*/
const char *Config::ConfigKeywords::IndividualCount		= "INDIVIDUAL_COUNT";
/**Number of families. This is actually the same as individual count	*/
const char *Config::ConfigKeywords::FamilyCount			= "FAMILY_COUNT";
const char *Config::ConfigKeywords::DatasetCount			= "DATASET_COUNT";
const char *Config::ConfigKeywords::AddDataset			= "DATASET";
const char *Config::ConfigKeywords::AddFamily				= "PED";
const char *Config::ConfigKeywords::PedContinuous			= "PED_CONT";
const char *Config::ConfigKeywords::AddRefPed				= "PED_REFERENCE";
const char *Config::ConfigKeywords::AddFamilyRandom		= "FLEXIPED";
const char *Config::ConfigKeywords::AddFamilyType			= "FAMTYPE";
const char *Config::ConfigKeywords::AddCC					= "CC";
const char *Config::ConfigKeywords::AddCont				= "CONT";
const char *Config::ConfigKeywords::AddContTails			= "CONT_TAILS";
const char *Config::ConfigKeywords::IncludeConfiguration 	= "INCLUDE";
/**
 * This is used to change the pedigree header from 10 to 6 columns
 */
const char *Config::ConfigKeywords::StandardPedigreeHeader = "USE_STD_PEDIGREE_HEADER";

/////Chromosom/block definition
const char *Config::ConfigKeywords::SetDefaultBlock 		= "DEFAULT_BLOCK";			///<Sets up the default block
const char *Config::ConfigKeywords::AddChromosome 		= "ADD_CHROMOSOME";			///<Start the definition of a new chromosome
const char *Config::ConfigKeywords::LoadChromosome		= "LOAD_CHROMOSOME";		///<Chromosome loci are defined in file
const char *Config::ConfigKeywords::LoadChromosomeXY		= "LOAD_CHROMOSOME_XY";		///<Chromosome loci defined in file (XY)
const char *Config::ConfigKeywords::AssocChromosome		= "ASSOC_CHROMOSOME";		///<Chromosome loci + association grid
const char *Config::ConfigKeywords::SeedChromosome		= "SEED_CHROMOSOME";
const char *Config::ConfigKeywords::AddBlock				= "ADD_BLOCK";				///<Add a block to the current chromosome

const char *Config::ConfigKeywords::FirstDropPoint 		= "FIRST_DROP_POINT";		///<Specify where you want to start dumping the gene pools
const char *Config::ConfigKeywords::DropFrequency 		= "DROP_FREQUENCY";			///<Specify how many generations are processed between the drops
const char *Config::ConfigKeywords::DropCount 			= "DROP_COUNT";				///<How many drops are made after the initial drop is made
const char *Config::ConfigKeywords::RandomizeAlleleFreq 	= "RANDOMIZE_ALLELE_FREQ";///<Turn on/off allele frequency randomization
const char *Config::ConfigKeywords::DoDumpGenerationZero 	= "DUMP_GENERATION_ZERO";

const char *Config::ConfigKeywords::TrueTypeFontName		= "FONT";				///<The font to be used
const char *Config::ConfigKeywords::BlockReportBuffer		= "LD_REPORT_BUFFER_SIZE";	///<The number of snps to be used on either side of a block in an LD plot
const char *Config::ConfigKeywords::CssFilename			= "CSS_FILENAME";		///<Override default stylesheet
const char *Config::ConfigKeywords::ClosePoolsBetweenAdvancements = "CLOSE_POOLS_BETWEEN_DROPS";
const char *Config::ConfigKeywords::UseOriginalCrossing	= "USE_ORIGINAL_CROSSING";	///<This is for debugging purposes only (and benchmarking)
const char *Config::ConfigKeywords::PedUseOriginalCross  = "PED_USE_ORIGINAL_CROSSING";
const char *Config::ConfigKeywords::PhasedPedigrees		= "PHASED_PEDIGREES";
const char *Config::ConfigKeywords::UseAltLD				= "USE_ALT_LD";

const char *Config::ConfigKeywords::DrawRSquared			= "DRAW_RSQUARED_PLOTS";
const char *Config::ConfigKeywords::DrawDPrime			= "DRAW_DPRIME_PLOTS";
const char *Config::ConfigKeywords::WriteLDReport			= "WRITE_LD_REPORT";

const char *Config::ConfigKeywords::MaxLdIndividualCount	= "FAST_LD_POOL_SIZE";
const char *Config::ConfigKeywords::PlotScanSize			= "FAST_LD_PLOT_SIZE";



const char *Config::ConfigKeywords::BinaryDatasets		= "BINARY_DATASETS";
const char *Config::ConfigKeywords::LocSelection			= "LOCUS_SELECTOR";
const char *Config::ConfigKeywords::AddRegion 			= "ADD_REGION";

const char *Config::ConfigKeywords::MaxThreadCount		= "MAX_THREAD_COUNT";
const char *Config::ConfigKeywords::OverideRecThreadCount = "OVERRIDE_REC_THREAD_COUNT";

const char *Config::ConfigKeywords::SimultaneousChroms	= "SIMULTANEOUS_CHROM";
const char *Config::ConfigKeywords::ThreadsPerChrom		= "THREADS_PER_CHROM";

const char *Config::ConfigKeywords::NoFastLD				= "NO_FAST_LD";
const char *Config::ConfigKeywords::TargetPopSize			= "TARGET_POP_SIZE";

const char *Config::ConfigKeywords::LociPerChromReported	= "MAX_LOCI_PER_CHROM_REPORTED";

const char *Config::ConfigKeywords::UseAdamEve 			= "INITIAL_POPULATION_ADAM_EVE";
const char *Config::ConfigKeywords::UseEden				= "INITIAL_POPULATION_EDEN";
const char *Config::ConfigKeywords::FounderCount			= "FOUNDER_COUNT";
const char *Config::ConfigKeywords::FounderDistortion		= "FOUNDER_DISTORTION";
const char *Config::ConfigKeywords::MaxRepeatCount		= "MAX_REPEAT_COUNT";
const char *Config::ConfigKeywords::ParentDistortion		= "PARENT_DISTORTION";
const char *Config::ConfigKeywords::ChildDistortion		= "CHILD_DISTORTION";

const char *Config::ConfigKeywords::WriteLOD				= "WRITE_LOD";



Config::Config() : currPool(NULL), doDumpGenerationZero(true), growthModelConfigured(false), configurationFilename("")  {
	growthModel = new LinearGrowth(100000, 0.0, 0.0);
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

void Config::SetRandomSeed( uint seed ) {
		generalSettings.seed = seed;
		Simulation::StatusModel::GAModel::seed = generalSettings.seed;
		Utility::Random::globalGenerator.Seed(generalSettings.seed);
		Simla::Random::setseed(seed);
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
		if (strcmp(value, "TIME") == 0) {
			struct timeval tv;
			struct timezone tz;
			gettimeofday(&tv, &tz);
			SetRandomSeed(tv.tv_sec*(tv.tv_usec<<0));
		}
		else {
			SetRandomSeed(atoi(value));
		}
	}
	else if (strcmp(key, ConfigKeywords::FamilyCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::FirstDropPoint) == 0)
		generalSettings.firstDropPoint = atoi(value);
	else if (strcmp(key, ConfigKeywords::DropFrequency) == 0)
		generalSettings.dropFrequency = atoi(value);
	else if (strcmp(key, ConfigKeywords::DropCount) == 0)
		generalSettings.dropCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::TargetPopSize) ==0)
		ChromPool::targetPop = atoi(value);
	else if (strcmp(key, ConfigKeywords::IndividualCount) == 0)
		datasetSettings.individualCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::DefineModel) ==0)
		statusSettings.DefineModel(line);
	else if (strcmp(key, ConfigKeywords::GenerateModels) == 0)
		GAModel::doRunGA = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::AddChromosome) == 0)
		AddChromosome(line);
	else if (strcmp(key, ConfigKeywords::LoadChromosome) == 0)
		LoadChromosome(line);
#ifdef USE_XY
	else if (strcmp(key, ConfigKeywords::LoadChromosomeXY) == 0)
		LoadChromosomeXY(line);
#endif //USE_XY
	else if (strcmp(key, ConfigKeywords::AssocChromosome) == 0)
		LoadAssocChromosome(line);
	else if (strcmp(key, ConfigKeywords::SeedChromosome) == 0)
		SeedChromosome(line);
	else if (strcmp(key, ConfigKeywords::AddBlock) == 0)
		AddBlock(line);
	else if (strcmp(key, ConfigKeywords::SetDefaultBlock) == 0)
		SetDefaultBlock(line);
	else if (strcmp(key, ConfigKeywords::RandomizeAlleleFreq) == 0)
		ChromPool::randomizeAlleleFreq = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::SimultaneousChroms) == 0)
		PoolManager::simultaneousChrom = atoi(value);
	else if (strcmp(key, ConfigKeywords::ThreadsPerChrom) == 0)
		PoolManager::threadsPerChrom = atoi(value);
	else if (strcmp(key, ConfigKeywords::NoFastLD) == 0)
		ChromPool::fastLD = false;
	else if (strcmp(key, ConfigKeywords::OffspringPerMating) == 0) {
		stringstream ss(line);
		string cmd;
		uint min, max;

		ss>>cmd>>min>>max;
		ChromPool::SetOffspringLimits( min, max );	
	}

	else if (strcmp(key, ConfigKeywords::MaxPoolSize) == 0) {
		if (growthModelConfigured) {
			cout<<"Attempting to set MAX_POOL_SIZE after the a model was defined (GROWTH_RATE). \nPlease make sure that the bounds are set prior to the growth rate.\n";
			abort();
		}
		GrowthRate::maxPoolSize = atoi(value);
	}
	else if (strcmp(key, ConfigKeywords::MinPoolSize) == 0) {
		if (growthModelConfigured) {
			cout<<"Attempting to set MIN_POOL_SIZE after the a model was defined (GROWTH_RATE). \nPlease make sure that the bounds are set prior to the growth rate.\n";
			abort();
		}

		GrowthRate::minPoolSize = atoi(value);
	}
	else if (strcmp(key, ConfigKeywords::GrowthRate) == 0)
		SetGrowthRate(line);
	else if (strcmp(key, ConfigKeywords::AddDataset) == 0)
		AddDataset(line);
	else if (strcmp(key, ConfigKeywords::MaxSnpsPerRow) == 0)
		Simulation::Visualization::ImageParameters::maxSnpsPerRow = atoi(value);
	else if (strcmp(key, ConfigKeywords::MaxSNPDistance) == 0)
		Simulation::Visualization::LdPlotter::maxSnpDistance = atoi(value);
	else if (strcmp(key, ConfigKeywords::GenerateOverviewLD) == 0)
		PoolManager::generateOverviewLD = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::GAFitnessThreshold) == 0) 
		Simulation::StatusModel::GAModel::fitnessThreshold = atof(value);
	else if (strcmp(key, ConfigKeywords::GAIterations) == 0)
		Simulation::StatusModel::GAModel::tries = atoi(value);
	else if (strcmp(key, ConfigKeywords::MappingFN) == 0)
		SetMappingFn(value);
	else if (strcmp(key, ConfigKeywords::LociPerChromReported) == 0)
		LocusSelection::maxReportedEntries = atoi(value);
	else if (strcmp(key, ConfigKeywords::StandardPedigreeHeader) ==0)
		Individual::StandardPedigreeHeader = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::DoDumpGenerationZero) == 0)
		doDumpGenerationZero = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::BlockReportSize) == 0) 
		ChromPool::numberOfBlocksToReport = atoi(value);
	else if (strcmp(key, ConfigKeywords::TrueTypeFontName) ==0) {
		stringstream s(line);
		string junk;
		s>>junk;
		Simulation::Visualization::ImageParameters::font = ParseFilename(s, "Font Filename");
	}
	else if (strcmp(key, ConfigKeywords::BlockReportBuffer) ==0)
		ChromPool::highBufferSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::CssFilename) == 0)
		ChromPool::cssFilename = value;
	else if (strcmp(key, ConfigKeywords::UseOriginalCrossing) == 0)
		ChromPool::UseOriginalCrossing = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::PedUseOriginalCross) == 0)
		PedigreeSample::UseOriginalCross = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::PhasedPedigrees) == 0)
		Individual::PhasedPedigrees = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::UseAltLD) == 0)
		ChromPool::UseAltLD = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::MaxLdIndividualCount) == 0)
		ChromPool::maxLDIndividuals = atoi(value);
	else if (strcmp(key, ConfigKeywords::PlotScanSize) == 0)
		ChromPool::plotScanSize = atoi(value);
	else if (strcmp(key, ConfigKeywords::DrawDPrime) == 0)
		ChromPool::writeDPrimePlots = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::DrawRSquared) == 0)
		ChromPool::writeRSquaredPlots = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::WriteLDReport) == 0)
		ChromPool::writeLdTextReport  = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::ClosePoolsBetweenAdvancements) == 0)
		generalSettings.closePoolsBetweenAdvancing = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::BinaryDatasets) == 0)
		datasetSettings.binaryDatasets = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::MinAlleleFreqThresh) == 0) {
		Individual::minAlFreqThreshold = atoi(value);
		if (Individual::minAlFreqThreshold > 0.0) {
			cout<<"!!!!!!!!!!!!!!!!!!!!!\n!! This is currently not recommended when generating datasets. Be warned\n";
		}
	}
	else if (strcmp(key, ConfigKeywords::IncludeConfiguration) == 0) {
		if (!LoadSettings(value))  {
			cout<<"Unable to load settings from the file: "<<value<<". Refusing to continue\n";
			abort();
		}	
	}
	else if (strcmp(key, ConfigKeywords::TemplateMatchAttempts) == 0)
		PedigreeTemplates::TemplatedPedigree::maxAttempts = atoi(value);
	else if (strcmp(key, ConfigKeywords::WriteLOD) == 0)
		PedigreeReferenceSample::DoWriteLOD = GetBoolean(value);
	else if (strcmp(key, ConfigKeywords::TemplateThresh) == 0)
		PedigreeTemplates::TemplatedPedigree::tolerance = atof(value);
	else if (strcmp(key, ConfigKeywords::LocSelection) == 0) 
		generalSettings.AddSelector( line );
	else if (strcmp(key, ConfigKeywords::AddRegion) == 0) 
		generalSettings.AddSelectorRegion( line );
	else if (strcmp(key, ConfigKeywords::UseAdamEve) == 0)
		ChromPool::UseAdamEve = true;
	else if (strcmp(key, ConfigKeywords::UseEden) == 0)
		ChromPool::UseEden = true;
	else if (strcmp(key, ConfigKeywords::FounderCount) == 0)
		ChromPool::FounderCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::FounderDistortion) == 0)
		ChromPool::FounderDistortion = atof(value);
	else if (strcmp(key, ConfigKeywords::MaxRepeatCount) == 0)
		ChromPool::MaxRepeatCount = atoi(value);
	else if (strcmp(key, ConfigKeywords::ParentDistortion) == 0)
		ChromPool::ParentDistortion = atof(value);
	else if (strcmp(key, ConfigKeywords::ChildDistortion) == 0)
		ChromPool::ChildDistortion = atof(value);
	else {
		cout<<"!! Unhandled line: "<<line<<"\n";
		success = false;
	}

	return success;
}

void Config::SetMappingFn(const char *val) {
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
	int locus=0, chromID=0;
	float al1=0.0, al2=0.0;

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
		cout<<"Correct syntax: \n\tALLELE_FREQ chromosome locus freq(A) freq(a)\n";
		abort();
	}
}

void Config::AddDataset(const char *line) {
	string cmd, type;
	uint columns = CountColumns(line);
	stringstream ss(line);

	//uint wordCount = CountColumns( line );

	ss>>cmd>>type;
	Sample *sample=NULL;

	static PedigreeSample* lastPedigree = NULL;			///<Help keep up with the last sample added

	if (strcmp(type.c_str(), ConfigKeywords::AddFamily)==0) {
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		string desc;
		ss>>desc>>genotypeError>>phenoError>>missingData;
		
		
		lastPedigree=new PedigreeSample(genotypeError, phenoError, missingData, desc.c_str());
		sample = lastPedigree;
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddFamilyType) == 0) {
		uint aff=0, uaff=0, rnd=0, count=0;
		ss>>aff>>uaff>>rnd>>count;
		if (lastPedigree) {
			lastPedigree->AddFamilyType(aff, uaff, count,  rnd);

			if (aff + uaff == 0 || count == 0) {
				cout<<"Invalid Family Type Definition:\n\t"<<line<<"\n";
				cout<<"Correct Syntax: \t\nDATASET PED affected_sibs unaffected_sibs extra_sibs genotype_error phenotype_error missing_data\n";
				abort();
			}
		}
		else {
			cout<<"Error: "<<cmd<<" "<<type<<" was used before a family was defined\n";
			abort();
		}

	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddFamilyRandom) == 0){ 
		uint maxFamSize=0;
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		string desc;

		ss>>desc>>maxFamSize>>genotypeError>>phenoError>>missingData;
		
		if (maxFamSize == 0) {
			cout<<"No upper limit to children count specified on FLEXIPED. \n\t"<<line<<"\n";
			cout<<"Correct Syntax: \t\nDATASET FLEXIPED affected_sibs unaffected_sibs genotype_error phenotype_error missing_data\n";
			abort();
		}
		
		sample=new PedigreeMixedSample(maxFamSize, genotypeError, phenoError, missingData, desc.c_str());
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddCC) == 0) {
		uint affCount = 0, unaffCount = 0;
		float  genotypeError=0.0, phenocopyError=0.0, missingData=0.0;
		string desc;

		if (columns != 8) {
			cout<<"Unexpected format: "<<line<<"\n";
			cout<<"-- Unable to configure case/control\n";	
	
			cout<<"Correct syntax: \n\tDATASET CC name [num aff] [num unaff] [genotype error] [phenocopy error] [missing data]\n";
			abort();
		}
		ss>>desc>>affCount>>unaffCount;
		ss>>genotypeError>>phenocopyError>>missingData;
	
		if (affCount < 1)
			cout<<"!! Due to the absence of affecteds present, no status column will be written\n";

		sample=new BasicSample(affCount, unaffCount, 0, genotypeError, phenocopyError, missingData, desc.c_str());
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddCont) == 0) {
		uint individuals = 0;
		float  genotypeError=0.0, phenocopyError=0.0, missingData=0.0;
		string desc;

		if (columns != 7) {
			cout<<"Unexpected format: "<<line<<"\n";
			cout<<"-- Unable to configure case/control\n";	
	
			cout<<"Correct syntax: \n\tDATASET CONT name [individual count] [genotype error] [phenocopy error] [missing data]\n";
			abort();
		}
		ss>>desc>>individuals;
		ss>>genotypeError>>phenocopyError>>missingData;
	
		if (individuals < 1) {
			cout<<"Correct syntax: \n\tDATASET CONT name [individual count] [genotype error] [phenocopy error] [missing data]\n";
			cout<<"-- Hint There must be at least 1 individidual in the dataset before drawing datasets.\n"<<line<<"\n";
			abort();
		}

		sample=new ContinuousSample(individuals, 0, 0, genotypeError, phenocopyError, missingData, desc.c_str());
	}	
	else if (strcmp(type.c_str(), ConfigKeywords::AddContTails) == 0) {
		uint lowerTail = 0, midSection = 0, upperTail = 0;
		float  genotypeError=0.0, phenocopyError=0.0, missingData=0.0;
		string desc;

		if (columns != 9) {
			cout<<"Unexpected format: "<<line<<"\n";
			cout<<"-- Unable to configure case/control\n";	
	
			cout<<"Correct syntax: \n\tDATASET CONT_TAILS name [num lower tail] [num mid] [num upper tail] [genotype error] [phenocopy error] [missing data]\n";
			abort();
		}
		ss>>desc>>lowerTail>>midSection>>upperTail;
		ss>>genotypeError>>phenocopyError>>missingData;
	
		if (lowerTail + midSection + upperTail < 1) {
			cout<<"Correct syntax: \n\tDATASET CONT_TAILS name [num lower tail] [num mid] [num upper tail] [genotype error] [phenocopy error] [missing data]\n";
			cout<<"-- Hint There must be at least 1 individidual in the dataset before drawing datasets.\n"<<line<<"\n";
			abort();
		}

		sample=new ContinuousSample(lowerTail, upperTail, midSection, genotypeError, phenocopyError, missingData, desc.c_str());
	}		
	else if (strcmp(type.c_str(), ConfigKeywords::PedContinuous) == 0) {
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		string desc;
		ss>>desc>>genotypeError>>phenoError>>missingData;
		
		lastPedigree=new PedigreeContinuousSample(genotypeError, phenoError, missingData, desc.c_str());
		sample = lastPedigree;
	}
	else if (strcmp(type.c_str(), ConfigKeywords::AddRefPed) == 0) {
		float genotypeError=0.0, phenoError=0.0, missingData=0.0;
		int replicates = 0;
		bool writeKinships = false;
		string kinshipsYesNo = "";
		string desc, source;
		
		ss>>source>>desc>>replicates>>kinshipsYesNo>>genotypeError>>phenoError>>missingData;
		writeKinships = GetBoolean(kinshipsYesNo.c_str());
		sample = new PedigreeReferenceSample(replicates, genotypeError, phenoError, missingData, desc.c_str());
		((PedigreeReferenceSample*)sample)->SetReference(source.c_str());
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
		cout<<"Correct syntax: \n\tDEFAULT_ALLELE_FREQ min_freq max_freq\n";

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
	growthModelConfigured = growthModel!=NULL;
/*	cout<<"Writing growth table for "<<cmd<<" growth\n";
	string filename = generalSettings.outputName + string(".") + growthModel->GetType() + ".csv";
	ofstream file(filename.c_str(), ios_base::out);
	growthModel->DiagramGrowth(file, 1, 5000, 1);
	*/
}


bool Config::LoadSettings(const char *filename) {
	configurationFilename = filename;
	char c  = '#';
	Utility::LineParser lp(c);
	return lp.Parse(filename, this) > 0;
}


//SEED_CHROMOSOME chr22.ped chr22.loc Chrom-22
void Config::SeedChromosome(const char *line) {
	string cmd, poolDataFilename, fmt, locFilename, label;
	int headerCols = 0;

	uint wordCount = CountColumns(line);

	if (wordCount < 5 || wordCount > 6) {
		cout<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		cout<<"Correct Syntax: \n\tSEED_CHROMOSOME PED/CC NumHeaderCols source_file loc_file [label]\n";
		cout<<"\tExample: SEED_CHROMOSOME PED 6 chr22.ped chr22.loc Chrom-22\n";
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

void Config::LoadAssocChromosome(const char *line) {
//	uint wordCount = CountColumns(line);
	string filename, cmd, label, gridFilename;
	stringstream ss(line);
	
	ss>>cmd;
	filename = ParseFilename(ss, "Chromosome File");
	gridFilename = ParseFilename(ss, "Grid File");
	
	ss>>label;
	while (!ss.eof()) {
		string s;
		ss>>s;
		label=label + " " + s;
	}

	if (filename.length() < 1) {
		stringstream ss;
		ss<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		ss<<"Correct syntax: \n\tLOAD_CHROMOSOME chromosome_file association_grid [label]\n";
		
		throw Utility::Exception::General(ss.str().c_str());
	}

	if (label.length() < 1)
		label = Utility::ExtractBaseFilename( filename.c_str());
	
	currPool = pools.AddChromosome(filename.c_str(), gridFilename.c_str(), label.c_str());
}
#ifdef USE_XY
void Config::LoadChromosomeXY(const char *line) {
	string filename, cmd, label;
	stringstream ss(line);
	
	ss>>cmd;
	filename = ParseFilename(ss, "XY Chromosome File");

	ss>>label;
	while (!ss.eof()) {
		string s;
		ss>>s;
		label=label + " " + s;
	}

	if (filename.length() < 1) {
		stringstream ss;
		ss<<"Misconfigured XY Chromosome: \n\t"<<line<<"\n";
		ss<<"Correct syntax: \n\tLOAD_CHROMOSOME_XY chromosome_file [label]\n";
		
		throw Utility::Exception::General(ss.str().c_str());
	}

	if (label.length() < 1)
		label = Utility::ExtractBaseFilename(filename.c_str());
	pools.AddChromosomeXY(filename.c_str(), label.c_str(), Utility::Random::globalGenerator);
}
#endif //USE_XY
void Config::LoadChromosome(const char *line) {
//	uint wordCount = CountColumns(line);
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
	label = StripTrailingWhitespace(label.c_str());
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

void Config::AddChromosome(const char *line) {
	uint blockCount;
	string cmd;
	string label;
	
	uint wordCount = CountColumns(line);

	if (wordCount < 2 || wordCount > 4) {
		cout<<"Misconfigured Chromosome: \n\t"<<line<<"\n";
		cout<<"Correct syntax: \n\tADD_CHROMOSOME block_count [label]\n";

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


void Config::AddBlock(const char *line) {
	if (currPool) {
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
		try {
			currPool->DefineBlock(min, max, blckMin, blckMax, snpMin, snpMax, prob);
		} catch (Utility::Exception::General& e) {
			cout<<"Error: "<<line<<"\n";
			cout<<e.GetErrorMessage()<<"\n";
			exit(1);
		}
	}
	else {
		cout<<"\nError. A Chromosome must be defined before block definitions are added to it\n";
		abort();
	}
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

bool Config::PostLoad() {
	cout<<"Intializing (this may take a few minutes)";
	uint estimatedGenerations = generalSettings.firstDropPoint + generalSettings.dropCount * generalSettings.dropFrequency;
	generalSettings.growthChartFilename = growthModel->DrawGrowthChart(generalSettings.outputName.c_str(), 1, estimatedGenerations, (uint)(estimatedGenerations *.01));
	bool success = pools.InitializeLoci(generalSettings.firstGeneration, generalSettings.outputName.c_str(), generalSettings.doLoad);
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
	cout<<".Done!\n";cout.flush();
	return success;
}

/**
 * @brief Check for missing data or impossible configuration settings
 */
bool Config::Validate() {
	bool isValid = true;

	

	return isValid;
}


PenetranceModel *Config::Status::LoadModel(const char *project, PoolManager *poolMgr) {
	poolMgr->PrepForSampling(project, poolMgr->GetCurrentGeneration(), model);

	model->Refresh(poolMgr);
	model->Load();
	return model;
}

PenetranceModel *Config::Status::InitModel(PoolManager *poolMgr, bool requireModel) {
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
			cout<<"There was no model specified. Please set up a disease model in the configuration file\n";
		else
			cout<<"Invalid model definition: "<<modelCfg<<". Does this file really exist?\n";
		abort();
	}
	return model;
	//model->Init(ss, poolMgr);
}

void Config::General::ReconcileSelections(PoolManager &pools) {
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
				cout<<"No locus was found to match, "<<label1<<". Please correct this and run again. \n"<<line.str()<<"\n";
				abort();
			}
			if (locusSelections.find(regionLbl.c_str()) == locusSelections.end() ) {
				cout<<"No matching selector was found, "<<regionLbl<<". Please correct this and run again. \n"<<line.str()<<"\n";
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

void Config::General::AddSelector(const char *selLine) {
	stringstream ss(selLine);

	string label;
	float mafTarget = 0.0, mafMin = 0.0, mafMax = 0.0;
	uint blTarget = 0, blMin = 0, blMax = 0;
	string command;	

	ss>>command>>label>>mafTarget>>mafMin>>mafMax>>blTarget>>blMin>>blMax;

	Simulation::Visualization::LocusSelection newSel(mafTarget, mafMin, mafMax, blTarget, blMin, blMax);
	if (locusSelections.find(label) != locusSelections.end()) {
		cout<<"Locus Selector label, "<<label<<" must be unique.";
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


}
