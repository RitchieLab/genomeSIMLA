//
// C++ Implementation: appconfiguration
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "appconfiguration.h"
#include "genomicregion.h"
#include <iomanip>
#include "pathway.h"
#include "magicdb.h"
#include "parislogs.h"
#include "utility/random.h"
#include "bin.h"
#include "knowledgebase.h"

using namespace std;

namespace Paris {

FileLogs FileLogs::logger;
bool FileLogs::WriteFeatureDetails			= false;
bool FileLogs::WriteDetailedFeatureReport	= false;

AppConfiguration::AppConfiguration()
{
}


AppConfiguration::~AppConfiguration()
{
}

void AppConfiguration::Init() {
	InitKey("VARIATION_FILENAME", "variations.bn",	"Variations data");
	InitKey("SETTINGS_DB", 			"bio-settings.cn","BioFilter data");
	InitKey("INCLUDE_KNOWLEDGE",	"2",					"List the various knowledge base (by KB ID) separated by spaces");
	InitKey("DATA_SOURCE",			"",					"filename containing snp,pvalue");
	InitKey("POPULATION",			"CEU",				"Set the population ID to match the population your data is drawn\n# from so that LD patterns can be used to expand the gene boundaries.");
	InitKey("REPORT_PREFIX", 		"",					"Prefix used for all reports");
	InitKey("REPORT_NAME",			"",					"Single word to describe data followed by optional long description\n# which can contain spaces (no returns, though). These will be used in some of the reports.");
	InitKey("LOAD_ALL_ALIASES",	"NO",					"Loads all aliases and generates a text report containing their associations");
	InitKey("HTML_REPORTS",			"NO",					"Write reports in html format (not all reports support HTML formatting");
	InitKey("BIN_SIZE",				"10000",				"Target number of features inside each bin. Paris will define bins to\n# be as close as possible to this number, but seldom will the count be exact.");
	InitKey("P_COUNT",				"1000",				"Number of permutations to be performed on each pathway.");
	InitKey("PATHWAY_SIG_THRESH",	"0.05",				"Threshold for determining the significance of a pathway (based on permutations)");
	InitKey("RESULTS_SIG_THRESH",	"0.05",				"Threshold for determining if a SNP is significant.");
	InitKey("GENE_BOUNDARY_EXTENSION", "50000",		"How many base pair locations up and down stream do we expand gene boundaries");
	InitKey("RANDOM_SEED",			"1371",				"Set the random seed used in permutations");
	InitKey("IGNORE_PVALUES_OF_ZERO", "ON",			"ON/OFF to ignore pvalues of zero. If they aren't ignored, they will be\n# counted as insignificant");
	InitKey("ALLOW_REDUNDANT_FEATURES", "OFF",		"ON/OFF to allow features common to multiple genes in the same pathway to\n# be counted multipe times");
	InitKey("COL_CHROMOSOME",		"1",					"Columnar location used for chromosome (1-22XY");
	InitKey("COL_RSID",				"2",					"Columnar location of the RS (can have rs prefix (caps or not) or just be a\n# numerical value");
	InitKey("COL_PVALUE",			"3",					"Columnar location of the pvalue to be used");
	InitKey("REFINEMENT_THRESHOLD_MIN", "0.03",		"The lower bound for borderline pvalues (set this to equal\n# REFINEMENT_THRESHOLD_MAX to not perform refinement)");
	InitKey("REFINEMENT_THRESHOLD_MAX", "0.07",		"The upper bound for borderline pvalues (set this to equal\n# REFINEMENT_THRESHOLD_MIN to not perform refinement)");
	InitKey("REFINEMENT_REP_COUNT",		"1000",		"The number of repeteated ptests performed when a pvalue is determined to be\n# borderline");
	InitKey("SHOW_ALL_ASSOCIATED_PATHWAYS", "OFF",		"When writing pathway investigation reports, do we show all pathways or only\n# the signficant ones?");
	InitKey("NEGATIVE_CONTROLS", "0",               "When set above 0, paris will perform that many negative control tests.\n");
}


void AppConfiguration::PrintSet(const char *key, vector<string>& settings) {
	vector<string>::iterator itr = settings.begin();
	vector<string>::iterator end = settings.end();

	cout<<setw(35)<<right<<key<<" : ";
	int count = 0;
	while (itr != end ) {
		if (count++>0)
			cout<<",";
		cout<<*itr++;
	}
	cout<<"\n";
}

void AppConfiguration::ReportConfiguration(std::ostream& os) {
	map<string, vector<string> >::iterator itr= strings.begin();
	map<string, vector<string> >::iterator end = strings.end();

	cout<<"-------------------- Configuration Parameters ----------\n";
	while (itr != end) {
		PrintSet(itr->first.c_str(), itr->second);
		itr++;
	}

	if (Analyzer::refinementThreshold.first == Analyzer::refinementThreshold.second)
		cout<<setw(35)<<right<<"Refinement"<<" : OFF\n";
}

void AppConfiguration::WriteConfiguration(std::ostream& os) {
	
}

void AppConfiguration::ExecuteConfiguration() {
	//Write all of our settings to the relevant variables in memory
	GenomicRegion::pvThreshold = GetDouble("RESULTS_SIG_THRESH");
	Utility::Random::globalGenerator.Seed(GetInteger("RANDOM_SEED"));
	Bin::rnd = Utility::Random::globalGenerator;
	//ParisResults::resultsDB = GetBoolean("RESULTS_DB");
	FileLogs::logger.Open (GetString("REPORT_PREFIX").c_str());
	Feature::IgnorePValueOfZero = GetBoolean("IGNORE_PVALUES_OF_ZERO");
	Gene::AllowRedundantFeatures = GetBoolean("ALLOW_REDUNDANT_FEATURES");
	Analyzer::refinementThreshold = std::pair<float, float>(GetDouble("REFINEMENT_THRESHOLD_MIN"), GetDouble("REFINEMENT_THRESHOLD_MAX"));
	KnowledgeBase::AnalysisRepCount = GetInteger("REFINEMENT_REP_COUNT");
}

}
