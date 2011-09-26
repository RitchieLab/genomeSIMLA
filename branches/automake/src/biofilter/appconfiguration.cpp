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
#include <iomanip>
#include "genegenemodel.h"
#include "bioapplication.h"
#include "kbgroup.h"
using namespace std;

namespace Biofilter {

AppConfiguration::AppConfiguration()
{
}


AppConfiguration::~AppConfiguration()
{
}

void AppConfiguration::Init() {
	InitKey("SETTINGS_DB",				"bio-settings.cn","BioFilter data");
	InitKey("MAX_GENE_COUNT",			"30",					"Max number of genes before we ignore the group");
	InitKey("SNPS_SOURCE",				"",					"The source file for the RS numbers in your dataset");
	InitKey("INCLUDE_GROUPS",			"",					"List the various groups (by group ID) separated by spaces");
	InitKey("INCLUDE_GROUP_FILE",		"",					"File containing group IDs to be the groups to be searched");
	//InitKey("MODEL_FILENAME", 		"NONE",				"Set the filename for the output model list (none writes to std-out)");
	InitKey("MODEL_BUFFER_INIT",		"1000000",			"Set the initial size of the model buffer. ");
	InitKey("MODEL_BUFFER_MAX",		"10000000",			"Set the upper limit to the buffer. Bigger -> faster, but must remain within\n# the limits of the hardware or could cause the application\n# to fail or become so slow that it will never complete.");
	InitKey("POPULATION",				"NO-LD",				"Set the population ID to match the population your data is drawn from so that\n# LD patterns can be used to expand the gene boundaries.");
	InitKey("GENE_BOUNDARY_EXTENSION", "0",			"How many base pair locations up and down stream do we expand gene boundaries (Only used if POPULATION is NO-LD)");
	InitKey("DISEASE_DEPENDENT",		"",					"Add one or more files containing disease dependent genes ");
	InitKey("PREFERRED_ALIAS",			"",					"User can specify aliases for genes (the alias must be present in the database");
	InitKey("REPORT_PREFIX", 			"",		 			"	Prefix used for all reports");
	InitKey("LOAD_ALL_ALIASES",		"NO",					"Loads all aliases and generates a text report containing their associations");
	InitKey("HTML_REPORTS",				"NO",					"Write reports in html format (not all reports support HTML formatting");
	InitKey("IMPLICATION_IDX_DUPLICATE_WEIGHT", "0.0",	"Weight applied to implication index for disease dependent groups are associated with both genes");
	InitKey("BINARY_MODEL_ARCHIVE",	"YES",				"Indicates whether to use a binary format instead of a text version");
	InitKey("DISEASE_DEPENDENT_LEVEL", "ALL_MODELS",	"ALL_MODELS, GROUP_LEVEL, DD_ONLY  These are used to determine selectivity of gene/gene models based on disease dependent relationships.");
	InitKey("COLLAPSE_ASSOCIATION_REPORT", "NO",			"When true, the associations reported cease where groups would generate models (if possible)");
	InitKey("ASSOCIATION_REPORT",		"NO",					"Produces association report");
	InitKey("ASSOCIATION_GRAPH",		"NO",					"Produces association graph input files");
	InitKey("SNP_REPORT",				"NO",					"Produces SNP report");
	InitKey("CLEANUP_RSIDS",			"ON",					"Attempt to identify merged and expired RS IDs in the SNP Source and modify them accordingly");

	//InitKey("EXPORT_SNP_MODELS",	"NO",					"Exports snp-snp models immediately after gene-gene model production");


	/**
	 */
}
void AppConfiguration::PrintSet(const char *key, vector<string>& settings, ostream& os) {
	vector<string>::iterator itr = settings.begin();
	vector<string>::iterator end = settings.end();

	os<<setw(35)<<right<<key<<" : ";
	int count = 0;
	while (itr != end ) {
		if (count++>0)
			os<<",";
		os<<*itr++;
	}
	os<<"\n";
}

void AppConfiguration::ReportConfiguration(std::ostream& os) {
	map<string, vector<string> >::iterator itr= strings.begin();
	map<string, vector<string> >::iterator end = strings.end();

	os<<"-------------------- Configuration Parameters ----------\n";
	while (itr != end) {
		PrintSet(itr->first.c_str(), itr->second, os);
		itr++;
	}	
}

void AppConfiguration::WriteConfiguration(std::ostream& os) {
}

void AppConfiguration::ExecuteConfiguration() {
	//Write all of our settings to the relevant variables in memory
	GeneGeneModel::DuplicateDD_Weight = GetDouble("IMPLICATION_IDX_DUPLICATE_WEIGHT");
	std::string setting = GetString("DISEASE_DEPENDENT_LEVEL");
	if (setting == "GROUP_LEVEL")
		Knowledge::KbGroup::DiseaseDependentRelationship = Knowledge::KbGroup::DD_GroupLevel;
	else if (setting == "DD_ONLY") {
		Knowledge::KbGroup::DiseaseDependentRelationship = Knowledge::KbGroup::DD_Only;
	}
	else
		Knowledge::KbGroup::DiseaseDependentRelationship = Knowledge::KbGroup::AllModels;

	Knowledge::KbGroup::CollapseAssociationReport = GetBoolean("COLLAPSE_ASSOCIATION_REPORT");
	BioApplication::geneExtension								= GetInteger("GENE_BOUNDARY_EXTENSION");
}

}
