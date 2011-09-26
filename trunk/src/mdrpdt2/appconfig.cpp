#include "appconfig.h"
#include "tstatistic.h"
#include "familyrepository.h"
#include "foldstatistic.h"
namespace MdrPDT {

ApplicationConfiguration *ConfigurationReader::appConfiguration = NULL;


void PedigreeConfiguration::Init() {
	InitKey("COMBO_START", 		"1",	"#Smallest order to be considered");
	InitKey("COMBO_STOP", 		"2", 	"#Largest order to be considered");
	InitKey("AFFECTED_VALUE",	"2", 	"#Value used to represent affected status");
	InitKey("UNAFFECTED_VALUE", "1",	"#Value used to represent the unaffected status");
	InitKey("MISSING_THRESHOLD","0.10",	"#Percentage of missing data at which point a locus will be ignored during analysis");

	InitKey("PTEST_COUNT",		"1000", "#Number of permutation tests to be run");
	InitKey("PTEST_SEED", 		"1397", "#Seed in shuffling data");
	InitKey("PTEST_SHORTCIRCUIT", "0", 	"#Short circuit ptests based on percentage (proves only that a model is NOT significant");
	InitKey("THREAD_COUNT",		"1",	"#Number of threads to be used during ptest evaluation");
	InitKey("INPUTFILE", 		"",	  	"#The pedigree dataset to be analyzed");
	InitKey("MERLIN_FORMAT",	"Yes",	"#Turn on merlin style (6 column header) If No, header is 10 columns");
	InitKey("MERLIN_DAT",		"",		"#Specify the location of dat file (merlin format only)");
	InitKey("REPORTMODELCOUNT", "5",	"#The number of models to be reported on");
	InitKey("REPORT_THRESHOLD", "0.0",  "#Threshold for reporting models\n");
	InitKey("EXPAND_TRIOS", 	"Yes", 	"#Expand trios into DSPS (otherwise they are ignored");
	InitKey("EXPAND_ALL",		"Yes",	"#Create 'virtual' sibling for each affected child"); 
	InitKey("PEDIGREE_EXCLUSIONS", "",	"#Indicate certain pedigrees to be ignored during analysis");
	InitKey("CROSSVALINTERVAL",	"5",	"#Number of cross validation folds the data is to be divided into. Use 1 for no cross validation");
	InitKey("VERBOSE_FOLDING",  "No",	"#Display the details of which families contribute to which fold");

	InitKey("MENDELIAN_ERROR_LEVEL", "1", "#Determine how to handle mendelian errors");
	InitKey("MENDELIAN_PEDIGREE_THRESHOLD", "0", "#Set the number of errors necessary to remove pedigree from analysis");
	InitKey("EXT_PEDIGREE", 	"pedigree", "# Extension to be appended onto the pedigree report");
	InitKey("EXT_DISTRIBUTION",	"dist", 	"# Extension to be appeneded for the distribution report");
	InitKey("EXT_CLEANED_DATA",	"ped",		"# Extension to be used for the cleaned data");
	InitKey("WRITE_CLEAN_DATAFILE",	"No",	"# Y/N write cleaned data to an accompanying file");
	InitKey("EXCLUDE_PEDIGREES", " ",	"# List of Pedigree IDs to be excluded from analysis");
	InitKey("EXCLUDE_LOCUS", 	" ", 	"# List of SNP numbers (index 1...N) indicating SNPs to be ignored from analysis");
}

void PedigreeConfiguration::ReportConfiguration(std::ostream& os) {
	int width=45;

	os<<setw(width)<<"Model Size: "<<Evaluation::EvaluationMethod::minModelSize<<"-"
			<<Evaluation::EvaluationMethod::maxModelSize<<"\n";
	int xvCount = GetInteger("CROSSVALINTERVAL");
	if (xvCount > 1)
		os<<setw(width)<<"Cross Validation Folds: "<<xvCount<<"\n";
	else
		os<<setw(width)<<"No Cross Validation: "<<"\n";

	os<<setw(width)<<"Permutation Tests: "<<GetInteger("PTEST_COUNT")<<"\n";
	os<<setw(width)<<"Random Number Seed: "<<GetInteger("PTEST_SEED")<<"\n";

}

void PedigreeConfiguration::ExecuteConfiguration() {
	Evaluation::EvaluationMethod::minModelSize 	= GetInteger("COMBO_START");
	Evaluation::EvaluationMethod::maxModelSize 	= GetInteger("COMBO_STOP");
	Individual::AffectedValue 					= GetChar("AFFECTED_VALUE");
	Individual::UnaffectedValue 				= GetChar("UNAFFECTED_VALUE");
	PedigreeRepository::verboseFoldingReport 	= GetBoolean("VERBOSE_FOLDING");

	Individual::UseMerlin						= GetBoolean("MERLIN_FORMAT");
	Evaluation::FoldEvaluation::reportThreshold				= GetDouble("REPORT_THRESHOLD");
}

void PedigreeConfiguration::WriteConfiguration(std::ostream& os) {
	int width = 35;
	os<<"########################### Input ###########################\n";
	os<<"# The dataset to be used\n";
	os<<setw(width)<<left<<"INPUTFILE"<<GetLine("INPUTFILE")<<"\n";
	os<<"# Input format. Merlin has 6 columns in ped file and an optional .dat file\n";
	os<<setw(width)<<left<<"MERLIN_FORMAT"<<GetLine("MERLIN_FORMAT")<<"\n";
	os<<"# Specify the location of the optional dat file (contains the labels for each locus\n";
	os<<setw(width)<<left<<"MERLIN_DAT"<<GetLine("MERLIN_DAT")<<"\n";
	os<<"# Set the value associated with affected status\n";
	os<<setw(width)<<left<<"AFFECTED_VALUE"<<GetLine("AFFECTED_VALUE")<<"\n";
	os<<"# Set the value associated with the unaffected status\n";
	os<<"# Individuals with neither Affected nor Unaffected will be ignored by the analysis\n";
	os<<setw(width)<<left<<"UNAFFECTED_VALUE"<<GetLine("UNAFFECTED_VALUE")<<"\n";
	os<<"# Ignore 0 or more pedigrees from the data\n";
	os<<setw(width)<<left<<"EXCLUDE_PEDIGREES"<<GetLine("EXCLUDE_PEDIGREES")<<"\n";
	os<<"# Ignore 0 or more SNPs (1...N) from the data\n";
	os<<setw(width)<<left<<"EXCLUDE_LOCUS"<<GetLine("EXCLUDE_LOCUS")<<"\n";
	os<<"# Maximum amount of missing data before a SNP is ignored from analysis\n";
	os<<setw(width)<<left<<"MISSING_THRESHOLD"<<GetDouble("MISSING_THRESHOLD")<<"\n";
	

	os<<"\n\n########################## Basic Settings ###################\n";
	os<<"# Describes the type of models of interest. Models consist of 1 or more SNPs.\n";
	os<<"# The minimum number of SNPs in a model to be investigate.\n";
	os<<"# Valid Settings: 1..COMBO_END\n";	
	os<<setw(width)<<left<<"COMBO_START"<<GetInteger("COMBO_START")<<"\n";
	os<<"# The maxmimum number of SNPs to be considered.\n";
	os<<"# Valid Settings: [COMBO_START..MAX_INT)\n";
	os<<setw(width)<<left<<"COMBO_STOP"<<GetInteger("COMBO_STOP")<<"\n";
	os<<"# Set the number of cross validation folds are used in analysis\n";
	os<<"# Recommend settings: 1, 5, 10 (1 is no cross validation)\n";
	os<<setw(width)<<left<<"CROSSVALINTERVAL"<<GetInteger("CROSSVALINTERVAL")<<"\n";
	os<<"# Maximum number of models to be reported\n";
	os<<setw(width)<<left<<"REPORTMODELCOUNT"<<GetInteger("REPORTMODELCOUNT")<<"\n";
	os<<"# Threshold for reporting models. At most REPORTMODELCOUNT will be reported (those are sorted according to the T-Statistic\n";
	os<<setw(width)<<left<<"REPORT_THRESHOLD"<<GetDouble("REPORT_THRESHOLD")<<"\n";
	os<<"# Write out details regarding the contents of each cross validation fold\n";
	os<<setw(width)<<left<<"VERBOSE_FOLDING"<<GetInteger("VERBOSE_FOLDING")<<"\n";
	os<<"\n\n";
	os<<"######################### Permutation Tests ##################\n";
	os<<"# Number of permutation runs to be executed. 1000 is recommended\n";
	os<<setw(width)<<left<<"PTEST_COUNT"<<GetInteger("PTEST_COUNT")<<"\n";
	os<<"# The seed associated with the tests. Each test gets a new seed\n";
	os<<setw(width)<<left<<"PTEST_SEED"<<GetInteger("PTEST_SEED")<<"\n";
	os<<"# Short circuit the permutation tests. See manual for explaination\n";
	os<<setw(width)<<left<<"PTEST_SHORTCIRCUIT"<<GetInteger("PTEST_SHORTCIRCUIT")<<"\n";
	os<<"# How many simultaneous threads will be run\n";
	os<<"# Each PTest can theoretically be run in it's own thread (Multiple threads\n";
	os<<"# will not benefit if no ptests are being performed)\n";
	os<<setw(width)<<left<<"THREAD_COUNT"<<GetInteger("THREAD_COUNT")<<"\n";

	os<<"\n\n";
	os<<"# Determine how to handle mendelian errors when encountered. \n";
	os<<"# Acceptable Values:\n";
	os<<"#   1 - Report errors, but do nothing\n";
	os<<"#   2 - Report errors, and zero out loci in families where genotyping error has been found\n";
	os<<"#   3 - Report errors and remove pedigrees where the number of genotyping errors exceeds threshold\n";
	os<<setw(width)<<left<<"MENDELIAN_ERROR_LEVEL"<<GetInteger("MENDELIAN_ERROR_LEVEL")<<"\n";
	os<<"# Set the threshold, if level is 3\n";
	os<<setw(width)<<left<<"MENDELIAN_PEDIGREE_THRESHOLD"<<GetInteger("MENDELIAN_PEDIGREE_THRESHOLD")<<"\n";

	os<<"\n\n\n";
	os<<"######################### Report Names ########################\n";
	os<<"# Pedigree report (genotypes and folding details)\n";
	os<<setw(width)<<left<<"EXT_PEDIGREE"<<GetLine("EXT_PEDIGREE")<<"\n";
	os<<"# Distribution report (each item from the distribution)\n";
	os<<setw(width)<<left<<"EXT_DISTRIBUTION"<<GetLine("EXT_DISTRIBUTION")<<"\n";
	os<<"# Write a copy of the 'cleaned' output to file\n";
	os<<setw(width)<<left<<"WRITE_CLEAN_DATAFILE"<<GetLine("WRITE_CLEAN_DATAFILE")<<"\n";
	os<<"# Extension of the file (above)\n";
	os<<setw(width)<<left<<"EXT_CLEANED_DATA"<<GetLine("EXT_CLEANED_DATA")<<"\n";
}


ApplicationConfiguration *ConfigurationReader::Load(const char *filename) {
	if (appConfiguration) {
		delete appConfiguration;
	}
	appConfiguration = new PedigreeConfiguration();
	appConfiguration->Init();
	if (filename) 
		appConfiguration->Parse(filename);
	appConfiguration->ExecuteConfiguration();
	if (filename)
		this->filename=filename;
	else 
		this->filename="";
	return appConfiguration;
}

ConfigurationReader::~ConfigurationReader() {
	if (appConfiguration)
		delete appConfiguration;
}


}


