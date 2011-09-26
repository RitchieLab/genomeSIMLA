//
// C++ Implementation: biofilter
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "biofilter.h"
#include <iostream>
#include "timestamp.h"
#include "ldcorrection.h"
#include <iomanip>
#include "genegenemodelreader.h"
#include "utility/strings.h"

using namespace std;
using namespace Utility;

namespace Biofilter {

Biofilter::Biofilter() : detailedCoverage(false), doLoadRegionAliases(false), doWriteModelCounts(false), ldConfiguration(""), silentRun(false) {
	action = BiofilterAction::NoAction;
}


Biofilter::~Biofilter()
{
	if (!silentRun)
		cerr<<bioApp.GetReportLog()<<reportLog.str();
}


void Biofilter::PrintBanner()  {
	if (!silentRun) {
		cerr<<"biofilter "<<APPMAJOR<<"."<<APPMINOR<<"."<<APPBUGFIX<<" ("<<BUILD_NUMBER<<") "<<BUILD_TYPE<<"  "<<BUILD_DATE<<"\n";
#ifdef USE_MPI
		cout<<"* This application is compiled to run on parallel computing systems using MPI\n";
##else
		cout<<"* (serial)\n";
#endif
		cerr<<"\nMarylyn Ritchie, William Bush and Eric Torstenson\nPlease forward any comments or errors to biofilter@chgr.mc.vanderbilt.edu\n\n";
	}
}

void Biofilter::PrintHelp() {
	silentRun = false;
	PrintBanner();
#ifdef USE_MPI
	cerr<<"usage: biofilter <configuration file> [ [command] ...] [ [parameter] ...]\n";
#else
	cerr<<"usage: biofilter <configuration file> \n";
#endif
	cerr<<"\nbiofilter is a standalone application for use in investigating possible SNP associations\n"
			"\tin a set of data which, through biological knowledge, might be worth investigating\n"; 
	cerr<<"Optional Commands Include:\n";
	cerr<<"\t-S [--sample-config]                       -- Print sample configuration to std-out\n";
	cerr<<"\t--report-gene-coverage gene-list-filename  -- Reports the snp count for the genes in genelist \n"
		 <<"\t                                              for the snps in snp-source\n";
	cerr<<"\t--filter-by-genes gene-list-filename       -- Lists gene name and rsid for each SNP inside each gene.\n"
		 <<"\t                                              gene names and rsids can both appear multiple times.\n";
	cerr<<"\t--inject-gene-information analysis-results chrom-col rs-col gene-list  \n"
		 <<"\t                                              Injects gene(s) at the end of the CSV file and writes the\n"
		 <<"\t                                              combined data to a new file.\n";
	cerr<<"\t--marker-info                              -- Reports each SNP and it's position/chromosome\n"
		 <<"\t                                              in a format acceptable by haploview\n";
	cerr<<"\t--model-report  model-list-filename        -- Generates a report containing the Genes and groups\n"
		 <<"\t                                              associated with each two snp model listed in the file\n";
	cerr<<"\t--list-genes                               -- Lists all genes that are covered by at least one SNP\n";
	cerr<<"\t--map-snps-to-gene gene-filename (or ALL)  -- Creates a comma separated file containing snp to gene mapping\n"
		 <<"\t                                              with the snp's relationship to the gene.\n";
	cerr<<"\t--snp-report                               -- \n";
	cerr<<"\nOptional Parameters Include:\n";
	cerr<<"\t-s [--snps] <snps filename>                -- Override the snp source file ont he commandline\n";
	cerr<<"\t-C [--coverage] <snps filename>            -- Add a file to coverage report list\n";
	cerr<<"\t-D [--detailed-coverage]                   -- (used with -C) adds extra details to coverage report\n";
	cerr<<"\t-X (--export-snp-models) min-impl max-models\n"
		 <<"\t                                           -- Writes Snp-Snp Models to file. This assumes a pre-existing \n"
		 <<"\t                                              gene-gene model file \n";
	cerr<<"\t-W [--write-models] min-impl max-models    -- Writes model list to files limitted to those with min-impl\n"
		 <<"\t                                              or greater with a target snp-snp model count of max-models\n";
	cerr<<"\t-m [--show-models]                         -- Writes contents of model file to screen in human\n"
	 	 <<"\t                                              readable form\n";
	cerr<<"\t-l [--load-ld] <ld filename>               -- Loads LD information from the file, filename, and\n"
		 <<"\t                                              adjusts the gene boundaries accordingly\n";
	cerr<<"\t-d [--disease-dependent] <filename>        -- Adds a meta group containing data from the file, filename\n";
	cerr<<"\t-G [--list-groups] [criteria]              -- Adds group search criteria and produces a list of\n"
		<<"\t                                               group IDs that match the criteria\n";
	cerr<<"\t-P [--list-populations]                    -- Lists all available Population based LD boundary options\n";
	cerr<<"\t-h [--html-reports] yes/no                 -- Turns HTML Reporting on/off\n";
	cerr<<"\t-b [--binary] yes/no                       -- Overrides binary setting in configuration file\n";
	cerr<<"\t-q [--quiet]                               -- Silences general output during processing. Reports and errors are still produced\n";
	cerr<<"\t-p [--set-population] pop                  -- Override the configurations population setting (NO-LD, CEUDP1.0, etc)\n";
	cerr<<"\t--optimize                                 -- Updates internal structures to allow faster access. This\n"
		 <<"\t                                              is usually done prior to release\n";
	cerr<<"\t--strip-optimization                       -- Strips the optimization out (this is helpful to allow data\n"
		 <<"\t                                              imports to run more quickly) \n";
	cerr<<"\t--ldspline ldconfig                        -- Imports LD-Spline variations using ldconfig as a guide\n";
	cerr<<"\t--fix-variations var-filename-path         -- Sets the path (and filename) to the appropriate variation file. ";
	cerr<<"\t                                              This should only be done if the file needs to be moved to a new location.\n";
}


int Biofilter::ParseCmd(int curr, int argc, char **argv) {
	int nextCmd = curr+1;
	if (strcmp(argv[curr], "--report-gene-coverage") == 0) 
		if (nextCmd < argc) {
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
			action = BiofilterAction::RunGeneCoverage;
		}
		else {
			action = BiofilterAction::ParseError;
			cerr<<"--report-gene-coverage must be followed the genelist filename\n";
			return -1;
		}
	else if (strcmp(argv[curr], "--map-snps-to-gene")== 0) {
		if (nextCmd < argc) {
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
			action = BiofilterAction::BiofilterMapSnpToGenes;
		}
		else {
			action = BiofilterAction::ParseError;
			cerr<<"--map-snps-to-gene must be followed by filename (contains SNP aliases on individual lines) or the keyword, ALL\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "--filter-by-genes")== 0) {
		if (nextCmd < argc) {
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
			action = BiofilterAction::FilterByGenes;
		} else {
			action = BiofilterAction::ParseError;
			cerr<<"--filter-by-genes must be followed by the filename containing genes (one on each line)\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "--inject-gene-information")==0) {
		if (nextCmd < argc +3) {
			cfg.SetValue("ANALYSIS_FILENAME", argv[nextCmd++]);
			cfg.SetValue("CHROM_COL", argv[nextCmd++]);
			cfg.SetValue("RS_COL", argv[nextCmd++]);
			cfg.SetValue("GENE_COVERAGE", argv[nextCmd++]);
			action = BiofilterAction::InjectGeneInformation;
		} else {
			action = BiofilterAction::ParseError;
			cerr<<"--inject-gene-information must be followed by the filename containing genes (one on each line)\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "--fix-variations") == 0) {
		if (nextCmd < argc) {
			action = BiofilterAction::SetVariationFilename;
			cfg.SetValue("VARIATIONS_FILENAME", argv[nextCmd++]);
		} else {
			action = BiofilterAction::ParseError;
			cerr<<"--fix-variations must be followed by the filename of the variations file to be used with this setting.\n";
		}
	}
	else if (strcmp(argv[curr], "--snp-report") == 0)
		cfg.SetValue("SNP_REPORT", "YES");

		//action = BiofilterAction::ProduceSnpReport;
	else if (strcmp(argv[curr], "--model-report") == 0) {
		action = BiofilterAction::RunModelReport;
		if (nextCmd < argc)
			cfg.SetValue("SNPS_SOURCE", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"--model-report must be followed by the model file\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "--list-associations") == 0)
		action = BiofilterAction::ListAssociations;
	else if (strcmp(argv[curr], "--graph-associations") == 0)
		action = BiofilterAction::GraphAssociations;
	else if (strcmp(argv[curr], "-D")==0 || strcmp(argv[curr], "--detailed-coverage")==0) 
		detailedCoverage = true;
	else if (strcmp(argv[curr], "-C")==0 || strcmp(argv[curr], "--coverage")==0) {
		if (nextCmd < argc)
			cfg.AppendValue("COVERAGE_SNPS", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-C (--coverage) must be followed by a snp filename\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-G")==0 || strcmp(argv[curr], "--list-groups")==0) {
		action = BiofilterAction::ListGroups;
		if (nextCmd < argc) {
			cfg.AppendValue("GROUP_SEARCH_CRITERIA", argv[nextCmd++]);
		}
	}
	else if (strcmp(argv[curr], "-P")==0 || strcmp(argv[curr], "--list-population-ids")==0) {
		action = BiofilterAction::ListPopulationIDs;
	}
	else if (strcmp(argv[curr], "-p")==0 || strcmp(argv[curr], "--set-population")==0) {
		if (nextCmd < argc)
			cfg.SetValue("POPULATION", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"--set-population must be followed by a population\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-d")==0 || strcmp(argv[curr], "--disease-dependent")==0) {
		if (nextCmd < argc)
			cfg.AppendValue("DISEASE_DEPENDENT", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"--disease-dependent must be followed by a filename\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-l")==0 || strcmp(argv[curr], "--load-ld")==0) {
		if (nextCmd < argc) {
			ldConfiguration = argv[nextCmd++];
		}
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-l (--load-ld) must be followed by a ld configuration filename\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-b")==0 || strcmp(argv[curr], "--binary")==0) {
		if (nextCmd < argc)
			cfg.SetValue("BINARY_MODEL_ARCHIVE", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-b (--binary) must be followed by an option: YES/NO\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-h")==0 || strcmp(argv[curr], "--html-reports")==0) {
		if (nextCmd < argc) 
			cfg.SetValue("HTML_REPORTS", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-h (--html-reports) must be followed by an option: YES/NO\n";
			return -1;
		}
	}		
	else if (strcmp(argv[curr], "-S")==0 || strcmp(argv[curr], "--sample-config")==0)
		action = BiofilterAction::PrintSampleConfig;
	else if (strcmp(argv[curr], "-p")==0 || strcmp(argv[curr], "--print-count-estimates")==0)
		doWriteModelCounts = true;
	else if (strcmp(argv[curr], "-x")==0 || strcmp(argv[curr], "--export-snp-models")==0) {
		if (nextCmd < argc-1) {
			cerr<<"EXPORT MODELS\n";
			action = BiofilterAction::ExportSnpModels;
			cfg.SetValue("MINIMUM_IMPLICATION_INDEX", argv[nextCmd++]);
			cfg.SetValue("MAX_SNP_MODEL_COUNT", argv[nextCmd++]);
		}
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-X (--export-snp-models) must be followed by 3 parameters: \n\t[gene-gene model filename] [minimum implication index] [max snp-snp model count]\n";
			return -1;
		}
	}
	else if (strcmp(argv[curr], "-W")==0 || strcmp(argv[curr], "--write-models")==0) {
		action = BiofilterAction::ProduceModels;
		if (curr < argc - 1 && argv[curr+1][0] != '-') {
			if (nextCmd < argc -1) {
				cfg.SetValue("EXPORT_SNP_MODELS", "YES");
				cfg.SetValue("MINIMUM_IMPLICATION_INDEX", argv[nextCmd++]);
				cfg.SetValue("MAX_SNP_MODEL_COUNT", argv[nextCmd++]);
			}
		}
	}
	else if (strcmp(argv[curr], "-m")==0 || strcmp(argv[curr], "--show-models")==0) {
		if (nextCmd < argc) {
			action = BiofilterAction::ListModels;
		}
	}
	else if (strcmp(argv[curr], "--marker-info")==0) 
		action = BiofilterAction::RunMarkerInfo;
	else if (strcmp(argv[curr], "--list-genes")==0)
		action = BiofilterAction::ListGenesSimple;
	else if (strcmp(argv[curr], "--strip-optimization") ==0)
		action = BiofilterAction::StripOptimization;
	else if (strcmp(argv[curr], "--optimize")==0) {
		action = BiofilterAction::Optimize;
	}
	else if (strcmp(argv[curr], "-q")==0 || strcmp(argv[curr], "--quite")==0) {
		silentRun = true;
	}
	else if (strcmp(argv[curr], "-s")==0 || strcmp(argv[curr], "--snps" )==0)
		if (nextCmd < argc)
			cfg.SetValue("SNPS_SOURCE", argv[nextCmd++]);
		else {
			action = BiofilterAction::ParseError;
			cerr<<"-s (snps) must be followed by the snps filename\n";
			return -1;
		}
	else if (strcmp(argv[curr], "--ldspline") == 0) {
		action = BiofilterAction::ImportLdSplines;
		ldConfiguration = argv[nextCmd++];
	}
	else {
		action = BiofilterAction::ParseError;
		cerr<<"Unknown argument: "<<argv[curr]<<"\n";
		return -1;
	}

	return nextCmd;
}

bool Biofilter::ParseCmdLine(int argc, char **argv) {

	//Test the DB connection
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
	if (argc < 2) {
		PrintHelp();
		return false;
	}
	int i=1;
	if (argv[1][0] != '-')
		LoadConfiguration(argv[i++]);
	//Work out any other cmd line arguments
	for (; i<argc && i>0;) {
		i=ParseCmd(i, argc, argv);
	}
	if (action == BiofilterAction::ParseError) {
		return false;
	}
	if (action == BiofilterAction::PrintSampleConfig) {
		cfg.Init();
		cfg.Write(cout);
		return false;
	}
	bioApp.SetReportPrefix(cfg.GetString("REPORT_PREFIX").c_str());
	bioApp.UseHtmlReports(cfg.GetBoolean("HTML_REPORTS"));
	bioApp.InitBiofilter(cfg.GetLine("SETTINGS_DB").c_str(), !silentRun);

	doLoadRegionAliases = cfg.GetBoolean("LOAD_ALL_ALIASES");

	if (!silentRun)
		cfg.ReportConfiguration(cerr);

	return true;
}

AppConfiguration *Biofilter::LoadConfiguration(const char *cfgFilename) {
	cfg.Init();
	cfg.SetValue("REPORT_PREFIX", ExtractBaseFilename(cfgFilename));
	if (cfgFilename) 
		cfg.Parse(cfgFilename);
	cfg.ExecuteConfiguration();
	if (cfgFilename)
		configFilename=cfgFilename;
	else 
		configFilename="";
	return &cfg;
}

vector<uint> Biofilter::LoadSNPs() {
	string snpFilename = cfg.GetLine("SNPS_SOURCE");
	set<uint> snps;
	vector<uint> snpList;
	if (snpFilename != "ALL") {
		ifstream file(snpFilename.c_str());
		if (!file.good()) {
			cerr<<"SNP data source, "<<snpFilename<<", appears unreadable. Unable to continue.\n";
			exit(1);
		}
		while (file.good()) {
			size_t rsID = 0;
			string id;
			file>>id;
			if (id.find("r") != string::npos || id.find("R") != string::npos)
				id=id.erase(0,2);
			rsID = atoi(id.c_str());
			if (rsID > 0) {
				snpList.push_back(rsID);
			}
		}
		if (!silentRun)
			cerr<<"\n"<<setw(35)<<snpFilename<<" : "<<snpList.size()<<" SNPs ";cout.flush();
	}
	vector<int> inclusions;
	string cleanReport			= "";
	if (cfg.GetBoolean("CLEANUP_RSIDS")) {
		bioApp.LoadConversion();
		cleanReport = GetReportFilename("snp-cleanup");
		//bioApp.CleanRSIDs(snpList, cleanReport.c_str(), !silentRun);
		//reportLog<<setw(45)<<right<<"SNP Cleanup Report: "<<cleanReport<<"\n";
	}
	snps.insert(snpList.begin(), snpList.end());
	//0 is used for missing SNPs. We don't want those
	snps.erase(0);
	//snpsRecorded can be larger than snps due to the fact that there might be more than a single SNP with the same RS number (we record them both)
	int snpsRecorded = bioApp.LoadVariations(snps, cleanReport.c_str(), !silentRun);
	if (!silentRun)
		cerr<<" ("<<snpsRecorded<<" matches in our database )\n";


	return snpList;
}

void Biofilter::InitGroupData() {
	vector<uint> inclusions;
	vector<string> groupInclusions;

	//Let's set up disease independant
	vector<string> diseaseIndependent;
	cfg.GetLines("DISEASE_DEPENDENT", diseaseIndependent);
	vector<string>::iterator lineItr 	= diseaseIndependent.begin();
	vector<string>::iterator linesEnd 	= diseaseIndependent.end();
	while (lineItr != linesEnd) 
		bioApp.AddUserDefinedGroup(lineItr++->c_str());

	cfg.GetLines("INCLUDE_GROUPS", groupInclusions);

	string groupFilename = cfg.GetString("INCLUDE_GROUP_FILE");
	if (groupFilename.length() > 0) {
		Utility::FileToArray conv;
		Utility::LineParser lp;
		lp.Parse(groupFilename.c_str(), &conv, false);
		groupInclusions.insert(groupInclusions.end(), conv.strings.begin(), conv.strings.end());
	}

	lineItr = groupInclusions.begin();
	linesEnd = groupInclusions.end();

	while (lineItr!=linesEnd) {
		stringstream ss(*lineItr);
		while (!ss.eof()) {
			string group = "";
			ss>>group;
			if (group!="") {
				inclusions.push_back(atoi(group.c_str()));
			}
		}
		lineItr++;
	}



	string regionAlias = cfg.GetLine("PREFERRED_ALIAS");
	string groupReport = GetReportFilename("dd-contents");
	ofstream gr(groupReport.c_str());
	if (regionAlias == "")
		regionAlias = "ALL";
	bioApp.LoadGroupData(cfg.GetInteger("MAX_GENE_COUNT"), inclusions, gr, cfg.GetLine("POPULATION").c_str(), regionAlias.c_str());
	bioApp.LoadRegionAliases(regionAlias);

	//bioApp.PrintSNPs(cout);
}

void Biofilter::SimpleGeneCoverage(){
	string db = cfg.GetLine("SETTINGS_DB");
	if (!silentRun)
		cerr<<"Loading All regions\n";

	bioApp.LoadRegions("ALL", cfg.GetLine("POPULATION").c_str());
	bioApp.LoadRegionAliases();
	if (!silentRun)
		cerr<<"Listing genes\n";
	bioApp.ListGenes(cout);
}
void Biofilter::InjectGeneInformation() {
	string db = cfg.GetLine("SETTINGS_DB");
	string analysisFile = cfg.GetLine("ANALYSIS_FILENAME");
	string genes = cfg.GetLine("GENE_COVERAGE");
	uint chrCol = cfg.GetInteger("CHROM_COL") - 1;
	uint rsCol = cfg.GetInteger("RS_COL") - 1;

	if (rsCol == (uint)-1) {
		cerr<<"Found zero based index in --inject-gene-information configuration. Adjusting both parameters to 1 based index.\n";
		chrCol+= 1;
		rsCol += 1;
	}

	vector<string> geneList;
	if (genes != "ALL")
		 geneList = FileToArray().OpenFile(genes.c_str());
	if (geneList.size() > 0 || genes == "ALL") {
		string analysisResults = GetReportFilename("-snp-analysis.csv");
		ofstream file(analysisResults.c_str());

		reportLog<<setw(45)<<right<<"SNP Analysis: "<<analysisResults<<"\n";
		string geneInc = "ALL";
		if (geneList.size() > 0)
			geneInc = "'" + Utility::ToString<vector<string> >(geneList, "','") + "'";

		bioApp.LoadRegions(geneInc, cfg.GetLine("POPULATION").c_str());
		SnpToGeneManager snpToGenes;
		bioApp.BuildSnpToGene(snpToGenes, geneList);
		ifstream anFile(analysisFile.c_str());
		char line[4096];
		vector<string> words;

		snpToGenes.PrintReport(cerr);

		while (!anFile.eof()) {
			anFile.getline(line, 4096);
			string strippedLine = string(line);
			size_t pos = strippedLine.rfind('\r');

			if (pos != string::npos)
				strippedLine.erase(pos);
			Utility::TokenizeString(strippedLine.c_str(), words, ",");
			uint rs = Utility::ExtractRsNumber(words[rsCol].c_str());
			char *chr = NULL;
			if (chrCol < (uint)-2)
				chr = (char*)words[chrCol].c_str();
			std::set<Knowledge::KbRegion*> genes = snpToGenes.GetGeneList(chr, rs);
			std::set<Knowledge::KbRegion*>::iterator itr = genes.begin();
			std::set<Knowledge::KbRegion*>::iterator end = genes.end();
			vector<string> geneHits;
			if (genes.size() == 0) {
				geneHits.push_back(" ");
			}

			while (itr != end) {
				geneHits.push_back((*itr)->CommonName() + ":" + (*itr)->Chromosome());
				//words.push_back(" http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + (*itr)->CommonName());
				itr++;
			}
			words.push_back(Utility::ToString<vector<string> >(geneHits, "|"));
			file<<Utility::ToString<vector<string> >(words, ",")<<"\n";
		}
	}
	else {
		cerr<<"The file, "<<genes<<", appears to be empty. No genes were found.\n";
	}
}

void Biofilter::SnpToGeneMapReport() {
	string db		= cfg.GetLine("SETTINGS_DB");
	string genes	= cfg.GetLine("GENE_COVERAGE");
	map<Knowledge::KbRegion::SnpType::Type, string> snpTypes;
	snpTypes[Knowledge::KbRegion::SnpType::Interior] = "Interior";
	snpTypes[Knowledge::KbRegion::SnpType::Exterior] = "Exterior";
	snpTypes[Knowledge::KbRegion::SnpType::Flanking] = "Flanking";
	snpTypes[Knowledge::KbRegion::SnpType::Unknown] = "Unknown";

	vector<string> genelist;

	if (genes != "ALL")
		genelist		= FileToArray().OpenFile(genes.c_str());
	if (genelist.size() > 0 || genes == "ALL") {
		string snpReport	= GetReportFilename("-snp-gene-map.csv");
		ofstream file(snpReport.c_str());
		reportLog<<setw(45)<<right<<"SNP / Gene Report: "<<snpReport<<"\n";
		string geneInc		= "'" + Utility::ToString<vector<string> >(genelist, "','") + "'";
		bioApp.LoadRegions(geneInc, cfg.GetLine("POPULATION").c_str());
		SnpToGeneManager snpToGenes;
		bioApp.BuildSnpToGene(snpToGenes, genelist);
		string chromosomes[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};

		for (int i=0; i<25; i++) {
			std::set<uint> snps = snpToGenes.GetSnpLists(chromosomes[i].c_str());
			std::set<uint>::iterator snpItr = snps.begin();
			std::set<uint>::iterator snpEnd = snps.end();

			while (snpItr != snpEnd) {
				std::set<Knowledge::KbRegion*> genes = snpToGenes.GetGeneList(chromosomes[i].c_str(), *snpItr);
				std::set<Knowledge::KbRegion*>::iterator regItr = genes.begin();
				std::set<Knowledge::KbRegion*>::iterator regEnd = genes.end();

				while (regItr != regEnd) {
					Knowledge::KbRegion* reg		= *regItr++;
					file<<chromosomes[i]<<",rs"<<*snpItr<<","<<reg->CommonName()<<",";
					Knowledge::KbRegion::SnpType::Type type = reg->GetSnpTypeByRS(*snpItr);
					file<<snpTypes[type]<<"\n";
				}
				snpItr++;
			}
		}
	}
}

void Biofilter::FilterByGenes() {
	string db = cfg.GetLine("SETTINGS_DB");
	string genes = cfg.GetLine("GENE_COVERAGE");

	vector<string> genelist;				// = FileToArray().OpenFile(genes.c_str());

	if (genes != "ALL")
		 genelist = FileToArray().OpenFile(genes.c_str());
	if (genelist.size() > 0 || genes == "ALL") {
		string snpReport = GetReportFilename("-snp-report.csv");
		ofstream file(snpReport.c_str());
		reportLog<<setw(45)<<right<<"SNP / Gene Report: "<<snpReport<<"\n";
		string geneInc = "'" + Utility::ToString<vector<string> >(genelist, "','") + "'";
		bioApp.LoadRegions(geneInc, cfg.GetLine("POPULATION").c_str());
		SnpToGeneManager snpToGenes;
		bioApp.BuildSnpToGene(snpToGenes, genelist);
		string chromosomes[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"};

		for (int i=0; i<25; i++) {
			std::set<uint> snps = snpToGenes.GetSnpLists(chromosomes[i].c_str());
			std::set<uint>::iterator snpItr = snps.begin();
			std::set<uint>::iterator snpEnd = snps.end();

			while (snpItr != snpEnd) {
				std::set<Knowledge::KbRegion*> genes = snpToGenes.GetGeneList(chromosomes[i].c_str(), *snpItr);
				file<<chromosomes[i]<<",rs"<<*snpItr;
				std::set<Knowledge::KbRegion*>::iterator regItr = genes.begin();
				std::set<Knowledge::KbRegion*>::iterator regEnd = genes.end();

				while (regItr != regEnd) {
					file<<","<<(*regItr)->CommonName();
					regItr++;
				}
				file<<"\n";
				snpItr++;
			}
		}



		//bioApp.ReportSNPsInGenes(genelist, file);
	}
	else {
		cerr<<"A problem was encountered opening file, "<<genes<<"\n";
		exit(1);
	}

}
void Biofilter::DetailGeneCoverage() {
	string db = cfg.GetLine("SETTINGS_DB");
	string genes = cfg.GetLine("GENE_COVERAGE");
	vector<string> snpFiles;
	if (!cfg.GetLines("COVERAGE_SNPS", snpFiles))
		cerr<<"Unable to find coverage files!\n";

	vector<string> genelist;
	if (genes != "ALL")
		 genelist = FileToArray().OpenFile(genes.c_str());
	if (genelist.size() > 0 || genes == "ALL") {
		string geneInc = "'" + Utility::ToString<vector<string> >(genelist, "', '") + "'";
		bioApp.LoadRegions(geneInc, cfg.GetLine("POPULATION").c_str());
		bioApp.DetailCoverage(genelist, snpFiles, detailedCoverage);
	}
	else {
		cerr<<"A problem was encountered opening file, "<<genes<<"\n";
		exit(1);
	}
}


void Biofilter::RunModelReport() {
/**
	Reporting::ModelReport *report;
	bool htmlReports = cfg.GetBoolean("HTML_REPORTS");
	if (htmlReports) {
		string filename = cfg.GetLine("REPORT_PREFIX") + "-model-report.html";
		report = new Reporting::ModelReportHTML(filename.c_str());
	}
	else
		report = new Reporting::ModelReport(cout);

	ifstream file(cfg.GetLine("SNPS_SOURCE").c_str());

	vector<uint> snps;
	while (file.good()) {
		string rs1, rs2;
		file>>rs1>>rs2;
		if (rs1.find("r") != string::npos || rs1.find("R") != string::npos) 
			rs1.erase(0,2);
		if (rs2.find("r") != string::npos || rs2.find("R") != string::npos) 
			rs2.erase(0,2);
		uint snp1 = atoi(rs1.c_str());
		uint snp2 = atoi(rs2.c_str());

		if (snp1 > 0 && snp2 > 0) {
			snps.push_back(snp1);
			snps.push_back(snp2);
			report->AddModel(snp1, snp2);
		}
	}

	//I'm worried that this is way too much logic for this function....
	//Maybe some of this should be done at the actual application level
	uint snpsLoaded = bioApp.InitSNPs(snps, cfg.GetLine("VARIATION_FILENAME").c_str());
	cerr<<"SNPS Loaded. "<<snps.size()<<" -> "<<snpsLoaded<<"\n";
	InitGroupData();
	bioApp.RunReport(report);
	delete report;
	 */
}
std::string Biofilter::GetReportFilename(const char *extension) {
	string prefix = cfg.GetLine("REPORT_PREFIX");
	string joint = ".";
	if (extension[0] == '-' || extension[0]=='.' || extension[0]=='_')
		joint = "";
	string filename = prefix + joint + string(extension);
	return filename;
}
void Biofilter::RunCommands() {
	bioApp.SetReportPrefix(cfg.GetLine("REPORT_PREFIX").c_str());
	switch (action) {
		case BiofilterAction::SetVariationFilename:
		{
			bioApp.SetVariationFilename(cfg.GetLine("VARIATIONS_FILENAME").c_str());
			return;
		}
		case BiofilterAction::PrintSampleConfig:
		{
			cfg.Write(cout);
			return;
		}
		case BiofilterAction::Optimize:
		{
			cerr<<"Optimizing\n";
			bioApp.PerformOptimization();
			return;
		}
		case BiofilterAction::StripOptimization:
			bioApp.StripOptimization();
			return;
		case BiofilterAction::ListGroups:
			{
				vector<string> keywords;
				cfg.GetLines("GROUP_SEARCH_CRITERIA", keywords);
				bioApp.ListGroupIDs(keywords);
			} return;
		case BiofilterAction::ListPopulationIDs:
			bioApp.ListPopulationIDs();
			return;
		case BiofilterAction::ListMetaGroups:
			{
				InitGroupData();
				bioApp.ListMetaGroups(cout);
			}
			return;
		case BiofilterAction::ListModels:
			{
				//InitGroupData();
				string filename = GetReportFilename(".gene-gene");
				string geneFilename = GetReportFilename(".genes");
				bool binaryArchive = cfg.GetBoolean("BINARY_MODEL_ARCHIVE");
				GeneGeneModelReader modelArchive(geneFilename.c_str(), filename.c_str(), binaryArchive);
				GeneGeneModelReader::iterator itr = modelArchive.begin();
				std::map<float, uint> counts;

				SnpModelCollection modelCollection;

				//This will consume ~2 gigabytes
				while (itr.GetModels(modelCollection, 10000000, 1) > 0) {
					SnpModelCollection::iterator sItr = modelCollection.begin();
					SnpModelCollection::iterator sEnd = modelCollection.end();

					while (sItr != sEnd) {
						SnpSnpModel* m = *sItr++;
						float score = m->ImplicationIndex();
						if (counts.find(score)==counts.end())
							counts[score] = 1;
						else
							counts[score]++;
						m->Write(cout, false);
						delete m;
					}
				}
				std::map<float, uint>::iterator citr = counts.begin();
				std::map<float, uint>::iterator cend = counts.end();
				cout<<"Model Generation Completed:\n"
					<<"Impl.\n"
					<<"Index\tCount\n";
				while (citr != cend) {
					cout<<setprecision(2)<<citr->first<<"\t"<<citr->second<<"\n";
					citr++;
				}
				//cerr<<"!!!!!!Skipping model report for now\n";
				//bioApp.ReportOnModels(filename.c_str());
			}
			return;
		case BiofilterAction::ExportSnpModels:
			{
				string filename = GetReportFilename("gene-gene");
				string geneFilename = GetReportFilename(".genes");
				string snpModelFilename = GetReportFilename(".snpsnp");
				bool binaryArchive = cfg.GetBoolean("BINARY_MODEL_ARCHIVE");
				uint minImplicationIndex = cfg.GetInteger("MINIMUM_IMPLICATION_INDEX");
				uint maxSnpModelCount = cfg.GetInteger("MAX_SNP_MODEL_COUNT");
				GeneGeneModelReader modelArchive(geneFilename.c_str(), filename.c_str(), binaryArchive);

				std::map<float, uint> counts = modelArchive.ArchiveSnpModels(snpModelFilename.c_str(), maxSnpModelCount, minImplicationIndex, binaryArchive);

				std::map<float, uint>::iterator citr = counts.begin();
				std::map<float, uint>::iterator cend = counts.end();
				cout<<"Model Generation Completed:\n"
					<<"Impl.\n"
					<<"Index\tCount\n";
				while (citr != cend) {
					cout<<setprecision(2)<<citr->first<<"\t"<<citr->second<<"\n";
					citr++;
				}
				reportLog<<setw(45)<<right<<"Snp Models: "<<snpModelFilename.c_str()<<"\n";

				//cerr<<"!!!!!!Skipping model report for now\n";
				//bioApp.ReportOnModels(filename.c_str());
				return;
			}
		case BiofilterAction::RunModelReport:
			RunModelReport();
			return;
		default:
			{}
	}

	//This is a special case, and shouldn't be done with anything else!
	if (ldConfiguration.length() > 0) {
		if (action == BiofilterAction::ImportLdSplines) 
			bioApp.ImportLdSplines(ldConfiguration.c_str());
		else
			bioApp.ImportLD(ldConfiguration.c_str(), cfg.GetLine("VARIATION_FILENAME").c_str());
		return;
	}
	
	//The rest of these needs this to be done before they begin
	vector<uint> snps = LoadSNPs();
	switch (action) {
		case BiofilterAction::RunMarkerInfo:
			bioApp.WriteMarkerInfo(cout, detailedCoverage);
			return;
		case BiofilterAction::RunGeneCoverage:
			DetailGeneCoverage();
			return;
		case BiofilterAction::ListGenesSimple:
			SimpleGeneCoverage();
			return;
		case BiofilterAction::BiofilterMapSnpToGenes:
			SnpToGeneMapReport();
			return;
		case BiofilterAction::FilterByGenes:
			FilterByGenes();
			return;
		case BiofilterAction::InjectGeneInformation:
			InjectGeneInformation();
			return;
		default:
			{}
	}
	
	InitGroupData();
	//if (doLoadRegionAliases)
	//	bioApp.LoadRegionAliases();
	int maxGeneCount = cfg.GetInteger("MAX_GENE_COUNT");

	if (doWriteModelCounts) {
		int maxGeneCount = cfg.GetInteger("MAX_GENE_COUNT");
		bioApp.SummarizeModelCounts(maxGeneCount);
	}
	if (cfg.GetBoolean("SNP_REPORT")) {
		bool writeHtml = cfg.GetBoolean("HTML_REPORTS");
		if (writeHtml) {
			string snpMissingFilename = GetReportFilename("_nogenes.txt");
			ofstream missing(snpMissingFilename.c_str());
			string snpReportFilename = GetReportFilename("_SNP_Report.html");
			ofstream file(snpReportFilename.c_str());
			bioApp.SnpReport(file, missing, snps, writeHtml);
			reportLog<<setw(45)<<right<<"SNP Report : "<<snpReportFilename<<"\n";
		}
		else {
			string snpReportFilename = GetReportFilename("_SNP_Report.txt");
			ofstream file(snpReportFilename.c_str());
			bioApp.SnpReport(cout, cout, snps, writeHtml);
		}
	}

	if (cfg.GetBoolean("ASSOCIATION_REPORT")) {
		bioApp.ListPresentAssociations(maxGeneCount);
	}

	if (cfg.GetBoolean("ASSOCIATION_GRAPH")) {
		bioApp.GraphPresentAssociations(GetReportFilename("dot").c_str(), maxGeneCount);
	}

	if (action == BiofilterAction::ProduceModels) {
		char tmpname[128];
		strcpy(tmpname, "modelsXXXXXX");
		int initBufferSize = cfg.GetInteger("MODEL_BUFFER_INIT");
		int maxBufferSize = cfg.GetInteger("MODEL_BUFFER_MAX");
		bool binaryArchive = cfg.GetBoolean("BINARY_MODEL_ARCHIVE");
		GeneGeneModelArchive repo(tmpname, initBufferSize, maxBufferSize, binaryArchive);

		//	ModelRepository repo(tmpname, initBufferSize, maxBufferSize);
		int maxGeneCount = cfg.GetInteger("MAX_GENE_COUNT");
		string geneGeneReport = GetReportFilename("gene-gene");
		string summaryReport  = GetReportFilename("-model-summary.txt");
		ofstream file(summaryReport.c_str());
		bioApp.ProduceModels(repo, file, maxGeneCount);
 		string geneFilename = GetReportFilename("genes");
		std::map<float, uint> counts = repo.Archive(geneFilename.c_str(), geneGeneReport.c_str());
		std::map<float, uint>::iterator itr = counts.begin();
		std::map<float, uint>::iterator end = counts.end();
		reportLog<<setw(45)<<right<<"Gene-Gene Model Summary: "<<summaryReport<<"\n";
		cout<<"Gene-Gene Model Summary (Snp-Snp Model Estimates)\n";
		cout<<setw(20)<<"Impl. Idx "<<setw(20)<<right<<"Count"<<"\n";
		cout<<setw(20)<<"-------------"<<setw(20)<<right<<"---------"<<"\n";
		while (itr != end) {
			cout<<setw(20)<<setprecision(2)<<itr->first<<setw(20)<<right<<itr->second<<"\n";
			itr++;
		}
		reportLog<<setw(45)<<right<<"Gene Mapping: "<<geneFilename<<"\n";
		reportLog<<setw(45)<<right<<"Gene-Gene Models: "<<geneGeneReport<<"\n";
		reportLog<<setw(45)<<right<<"Summary Report: "<<summaryReport<<"\n";
	}
	if (cfg.GetBoolean("EXPORT_SNP_MODELS")) {
		string filename = GetReportFilename("gene-gene");
		string geneFilename = GetReportFilename("genes");
		string snpModelFilename = GetReportFilename("snpsnp");
		bool binaryArchive = cfg.GetBoolean("BINARY_MODEL_ARCHIVE");
		uint minImplicationIndex = cfg.GetInteger("MINIMUM_IMPLICATION_INDEX");
		uint maxSnpModelCount = cfg.GetInteger("MAX_SNP_MODEL_COUNT");
		GeneGeneModelReader modelArchive(geneFilename.c_str(), filename.c_str(), binaryArchive);

		std::map<float, uint> counts = modelArchive.ArchiveSnpModels(snpModelFilename.c_str(), maxSnpModelCount, minImplicationIndex, binaryArchive);

		std::map<float, uint>::iterator citr = counts.begin();
		std::map<float, uint>::iterator cend = counts.end();
		cout<<"\nSnp-Snp Model Generation Summary:\n"
			<<setw(20)<<"Impl."<<"\n"
			<<setw(20)<<right<<"Index "<<setw(20)<<"Count"<<"\n";
		cout<<setw(20)<<"-------------"<<setw(20)<<right<<"---------"<<"\n";
		while (citr != cend) {
			cout<<setw(20)<<setprecision(2)<<citr->first<<setw(20)<<citr->second<<"\n";
			citr++;
		}
		reportLog<<setw(45)<<right<<"Snp Models: "<<snpModelFilename.c_str()<<"\n";

	}
	cout<<"\n";

}



string Biofilter::GetReportPrefix() {
	string prefix = cfg.GetLine("REPORT_PREFIX");
	if (prefix == "")
		prefix = configFilename;
	return prefix;
}

}
int main(int argc, char *argv[])	{
	string cfgFilename;

	Biofilter::Biofilter app;					///<The application object
	if (!app.ParseCmdLine(argc, argv))
		exit(1);
	//Performs any commands
	app.RunCommands();
	
  	return EXIT_SUCCESS;
}
