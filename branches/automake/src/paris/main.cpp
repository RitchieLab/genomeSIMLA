/* 
 * File:   main.cpp
 * Author: torstees
 * 
 * Created on January 6, 2010, 12:54 PM
 */

#include "main.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "timestamp.h"
#include "utility/strings.h"
#include "parislogs.h"
using namespace std;
using namespace Utility;

namespace Paris {


Main::Main() : action(ParisAction::NoAction) {
}

Main::Main(const Main& orig) : action(orig.action) {
}

Main::~Main() {
}


void Main::PrintBanner()  {
	cerr<<"          *******************************************************\n";
	stringstream ss;
	ss<<APPMAJOR<<"."<<APPMINOR<<"."<<APPBUGFIX<<"b ("<<BUILD_NUMBER<<") "<<BUILD_TYPE;
	cerr<<"          * paris        "<<setw(35)<<right<<ss.str()<<"    *\n"
		 <<"          *                                                     *\n"
	    <<"          *              "<<setw(35)<<right<<BUILD_DATE<<"    *\n";
#ifdef USE_MPI
	cerr<<"* This application is compiled to run on parallel computing systems using MPI\n";
##else
	cerr<<"* (serial)\n";
#endif
	cerr<<"          * Brian Yaspan, Jonathan Haines, Marylyn Ritchie and  *\n"
		 <<"          * Eric Torstenson                                     *\n"
	    <<"          * Comments/errors to TBD@chgr.mc.vanderbilt.edu       *\n"
		 <<"          *******************************************************\n";
}

void Main::PrintHelp() {
	PrintBanner();
#ifdef USE_MPI
	cerr<<"usage: paris <configuration file> [ [command] ...] [ [parameter] ...]\n";
#else
	cerr<<"usage: paris <configuration file> [options]\n";
#endif
	cerr<<"Optional Commands Include:\n";
	cerr<<"  -S [--sample-config]                 -- Print sample configuration to std-out\n";
	cerr<<"  -I [--investigate-pathways] filename -- Generate gene reports for each gene found\n"
		 <<"                                          in the list of groups specified by filename\n";
	cerr<<"  -K [--list-knowledge]                -- Lists the knowledge bases with their keys\n"
		 <<"                                          (keys are used in the 'INCLUDE_GROUPS'\n"
		 <<"                                          configuration parameter)\n";
	cerr<<"  -B [--bin-report]                    -- Create bin report\n";
	cerr<<"\nOptional Parameters Include:\n";
	cerr<<"  -P [--list-populations]              -- Lists all available Population based LD\n"
		 <<"                                          boundary options\n";
}


StringArray Main::LoadGroupFile(const char *groups) {
	Utility::LineParser lp;
	Utility::FileToArray fta;
	lp.Parse(groups, &fta);
	return fta.strings;
}

int Main::ParseCmd(int curr, int argc, char **argv) {
	int nextCmd = curr+1;
/*	if (strcmp(argv[curr], "-G")==0 || strcmp(argv[curr], "--list-groups")==0) {
		action = ParisAction::ListGroups;
		if (nextCmd < argc) {
			cfg.AppendValue("GROUP_SEARCH_CRITERIA", argv[nextCmd++]);
		}
	}
*/	if (strcmp(argv[curr], "-K")==0 || strcmp(argv[curr], "--list-knowledge")==0) {
		action = ParisAction::ListGroups;
	}
	else if (strcmp(argv[curr], "-I")==0 || strcmp(argv[curr], "--investigate-pathways")==0) {
		action = ParisAction::InvestigatePathways;
		if (nextCmd < argc) {
			cfg.AppendValue("PATHWAY_INVESTIGATION_FILENAME", argv[nextCmd++]);
		}
		else {
			action = ParisAction::ParseError;
			cerr<<"-I (--investigate-groups) must be followed by the filename which holds the groups to be investigated\n";
		}
	}
	else if (strcmp(argv[curr], "-P")==0 || strcmp(argv[curr], "--list-population-ids")==0) {
		action = ParisAction::ListPopulationIDs;
	}

	else if (strcmp(argv[curr], "-S")==0 || strcmp(argv[curr], "--sample-config")==0)
		action = ParisAction::PrintSampleConfig;
	else if (strcmp(argv[curr], "-B")==0 || strcmp(argv[curr], "--bin-report")==0)
		action = ParisAction::WriteBinReport;
	else {
		action = ParisAction::ParseError;
		cerr<<"Unknown argument: "<<argv[curr]<<"\n";
		return -1;
	}

	return nextCmd;
}

bool Main::ParseCmdLine(int argc, char **argv) {

	//Test the DB connection
#ifdef USE_MPI
	MPI::Init(argc, argv);
#endif
	if (argc < 2) {
		PrintHelp();
		return false;
	}
	PrintBanner();
	int i=1;
	if (argv[1][0] != '-')
		LoadConfiguration(argv[i++]);
	//Work out any other cmd line arguments
	for (; i<argc && i>0;) {
		i=ParseCmd(i, argc, argv);
	}
	if (action == ParisAction::ParseError) {
		return false;
	}
	if (action == ParisAction::PrintSampleConfig) {
		cfg.Init();
		cfg.Write(cout);
		return false;
	}
	parisApp.SetReportPrefix(cfg.GetString("REPORT_PREFIX").c_str());
	parisApp.ReportName(cfg.GetString("REPORT_NAME").c_str());

	parisApp.InitKnowledge(cfg.GetLine("SETTINGS_DB").c_str());

	//doLoadRegionAliases = cfg.GetBoolean("LOAD_ALL_ALIASES");

	cfg.ReportConfiguration(cout);

	return true;
}

AppConfiguration *Main::LoadConfiguration(const char *cfgFilename) {
	cfg.Init();
	cfg.SetValue("REPORT_PREFIX", Utility::ExtractBaseFilename(cfgFilename));
	if (cfgFilename)
		cfg.Parse(cfgFilename);
	cfg.ExecuteConfiguration();
	if (cfgFilename)
		configFilename=cfgFilename;
	else
		configFilename="";
	return &cfg;
}


string Main::GetReportPrefix() {
	string prefix = cfg.GetLine("REPORT_PREFIX");
	if (prefix == "")
		prefix = configFilename;
	return prefix;
}


void Main::RunCommands() {
	parisApp.SetReportPrefix(cfg.GetLine("REPORT_PREFIX").c_str());
	switch (action) {
		case ParisAction::PrintSampleConfig:
		{
			cfg.Write(cout);
			return;
		}
		case ParisAction::ListGroups:
			{
				parisApp.ListGroupIDs();
//				parisApp.ListGroupIDs(keywords);
			} return;


/*
		case ParisAction::ListPopulationIDs:
			parisApp.ListPopulationIDs();
			return;
		case ParisAction::ListMetaGroups:
			{
				InitGroupData();
				parisApp.ListMetaGroups(cout);
			}
			return;
		 */
		default: {}
	}

	multimap<string, uint> snps;
	vector<ParisApp::DataPoint> data;
	LoadDataPoints(data, snps);
	string ldPopulation					= cfg.GetLine("POPULATION");
	uint geneExpansion					= cfg.GetInteger("GENE_BOUNDARY_EXTENSION");
	StringArray knowledgeBases;
	cfg.GetLines("INCLUDE_KNOWLEDGE", knowledgeBases);
	parisApp.InitKB(ldPopulation.c_str(), geneExpansion, knowledgeBases);

	uint binSize							= cfg.GetInteger("BIN_SIZE");
	uint pCount								= cfg.GetInteger("P_COUNT");
	float pathSig							= cfg.GetDouble("PATHWAY_SIG_THRESH");
	float dataSig							= cfg.GetDouble("RESULTS_SIG_THRESH");
	parisApp.InitData(data, binSize);
	std::string vitalsReport			= GetReportPrefix() + "-kb-stats.txt";
	parisApp.WriteKbVitalsReport(vitalsReport.c_str());
	cerr<<"KB Statistics Report: "<<vitalsReport<<"\n";

	if (action == ParisAction::InvestigatePathways) {
		StringArray groups				= LoadGroupFile(cfg.GetString("PATHWAY_INVESTIGATION_FILENAME").c_str());
		StringArray::iterator itr		= groups.begin();
		StringArray::iterator end		= groups.end();

		while (itr != end) {
			parisApp.InvestigatePathway(pCount, dataSig, pathSig, itr++->c_str(), cfg.GetBoolean("SHOW_ALL_ASSOCIATED_PATHWAYS"));
		}
	}
	else if (action == ParisAction::WriteBinReport) {
		std::string filename = GetReportPrefix() + "-bin-report.txt";
		parisApp.WriteBinReport(filename.c_str());
		cerr<<"Bin Report: "<<filename<<"\n";
	}
	else {
		uint ncCount = cfg.GetInteger("NEGATIVE_CONTROLS");
		if (ncCount > 0) {
			std::map<std::string, vector<uint> > sigPathways;
			for (int i=0; i<ncCount; i++) {
				FileLogs::logger.OpenFinalReport(i+1);
				std::map<std::string, uint> results = parisApp.AnalyzeNegativeControl(pCount, dataSig, pathSig);
				std::map<std::string, uint>::iterator itr = results.begin();
				std::map<std::string, uint>::iterator end = results.end();
				while (itr != end) {
					cout<<"\t"<<itr->first<<"\t"<<itr->second<<"\n";
					sigPathways[itr->first].push_back(itr->second);
					itr++;
				}
			}

			std::map<std::string, vector<uint> >::iterator kb = sigPathways.begin();
			std::map<std::string, vector<uint> >::iterator kbEnd = sigPathways.end();
			float total=0.0;
			int count = 0;
			cout<<"Final Tally of the Negative Controls:\n";
			while (kb!=kbEnd) {
				cout<<"--Knowledge Base: "<<kb->first<<"\n";
				std::vector<uint>::iterator itr = kb->second.begin();
				std::vector<uint>::iterator end = kb->second.end();

				while (itr != end) {
					cout<<"\t"<<*itr<<"\n";
					total += (float)*itr;
					count++;
					itr++;
				}
				kb++;
				cout<<"------------\nAvg:\t"<<(total/(count))<<"\n";
			}
		} else {
			FileLogs::logger.OpenFinalReport();
			parisApp.RunAnalysis(pCount, dataSig, pathSig);
		}
	}

}


void Main::LoadDataPoints(vector<ParisApp::DataPoint>& data, std::multimap<string, uint >& snps) {
	string filename	= cfg.GetLine("DATA_SOURCE");		///< chr, rs, pos
	uint chrColumn		= cfg.GetInteger("COL_CHROMOSOME") - 1;	///< 0
	uint rsColumn		= cfg.GetInteger("COL_RSID") - 1;			///< 1
	uint pvalColumn	= cfg.GetInteger("COL_PVALUE") - 1;		///< 3



	if (filename == "") {
		cerr<<"There is no data source selected. Paris requires a pvalues to process pathway information.\n";
		exit(1);
	}
	char line[4096];
	ifstream file(filename.c_str());

	if (!file.good()) {
		cerr<<"Unable to open source data file, "<<filename<<". Aborting\n";
		exit(1);
	}
	uint lineNumber	= 0;
	while (file.good()) {
		file.getline(line, 4096);
		lineNumber++;
		vector<string> tokens;
		if (Utility::TokenizeString(line, tokens, "\t, ") > 2) {
			string chr	= tokens[chrColumn];
			if (chr == "23")
				chr = "X";
			else if (chr == "24")
				chr = "Y";
	
			string rs	= tokens[rsColumn];
			//string p		= tokens[2];
			string pvalue = tokens[pvalColumn];

			if (rs.find("r") != string::npos || rs.find("R") != string::npos)
					rs=rs.erase(0,2);

			if (rs[0] >= '0' and rs[0] <= '9') {
				uint rsID = atoi(rs.c_str());
				if ((chr == "X" || chr == "Y" || (atoi(chr.c_str()) > 0 && atoi(chr.c_str()) < 23))) {
	//	EST-RS>0			if (rs > 0.00) {
						data.push_back(ParisApp::DataPoint(atoi(rs.c_str()), chr.c_str(), atof(pvalue.c_str()) ));
						snps.insert(pair<string, uint>(chr, rsID));
	//				}
				}
				else {
					cerr<<"Ignoring line #"<<lineNumber<<". Chromosome "<<chr<<" is not recognized.\n";
				}
			}
			else {
				cerr<<"Ignoring SNP #"<<rs<<". Paris is currently only able to recognize RS numbers.\n";
			}
		}
	}
	cout<<"\n"<<setw(35)<<filename<<" : "<<snps.size()<<" SNPs ";cout.flush();

	parisApp.InitSNPs(snps);
	cout<<" ("<<snps.size()<<" matches in our database )\n";

}


}
#ifndef TEST_APP
int main(int argc, char **argv) {
	string cfgFilename;

	Paris::Main app;					///<The application object
	if (!app.ParseCmdLine(argc, argv))
		exit(1);
	//Performs any commands
	app.RunCommands();

  	return EXIT_SUCCESS;
}
#endif //TEST_APP
