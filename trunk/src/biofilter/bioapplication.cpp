//
// C++ Implementation: bioapplication
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "bioapplication.h"
#include "utility/strings.h"
#include "utility/exception.h"
#include "kbgroup.h"
#include "ldcorrection.h"
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <iomanip>
#include "modelreport.h"
#include <iomanip>
#include "utility/strings.h"
#include "ldsplineimporter.h"

using namespace soci;
using namespace std;
using namespace Spline;


namespace Biofilter {

using namespace Knowledge;

bool BioApplication::CrossGroupModelGen	= false;
uint BioApplication::geneExtension			= 0;
bool BioApplication::PurgeOldIDs				= false;
BioApplication::BioApplication(const char *prefix, bool htmlReports) 
		: reportPrefix(prefix), htmlReports(htmlReports), populationDesc("") { }


BioApplication::~BioApplication()
{
	map<uint, KbMetaGroup*>::iterator mitr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator mend = metagroups.end();
	while (mitr != mend)  
		delete (mitr++)->second;

	GroupLookup::iterator gItr = groups.begin();
	GroupLookup::iterator gEnd = groups.end();

	while (gItr != gEnd) 
		delete (gItr++)->second;


}
void BioApplication::SetReportPrefix(const char *prefix) {
	reportPrefix = prefix;
}

void BioApplication::UseHtmlReports(bool doUse) {
	htmlReports = doUse;
}

KbRegion *BioApplication::GetRegion(uint geneID) {
	KbRegion *region = NULL;
	
	if (regions.find(geneID) != regions.end())
		region = regions[geneID];

	return region;
}
string BioApplication::GetReportLog() {
	return reportLog.str();
}

void BioApplication::ImportLdSplines(const char *ldConfiguration) {
	LdSplineImporter splineMgr;
	splineMgr.LoadConfiguration(ldConfiguration);
	splineMgr.Process(sociDB);
}

void BioApplication::ImportLD(const char *ldConfiguration, const char *variationFilename) {
	LdCorrection ldImport;
	ldImport.LoadConfiguration(sociDB, ldConfiguration);
	ldImport.Process(sociDB, variationFilename);
}


void BioApplication::GraphPresentAssociations(const char *filename, uint maxGeneCount) {
	map<uint, KbMetaGroup*>::iterator itr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator end = metagroups.end();

	while (itr != end) {
		if (itr->second->GetGroupCount() > 0 && !itr->second->IsDiseaseDependent()) {
			string filename = reportPrefix + string("-") + itr->second->CommonName() + ".dot";
			ofstream file(filename.c_str());
			file<<"digraph G{\n";
			  //<<"\tsize=\"7,10\"\n"
			  //<<"\tpage=\"8.5,11\"\n"
			  //<<"\tcenter=\"\"\n"
			  file<<"\tnode[colorscheme=set312,style=\"filled\"]\n";
			itr->second->GraphAssociations(file, maxGeneCount);
			file<<"}";
			reportLog<<setw(45)<<right<<"Graph Association: "<<filename<<"\n";
		}
		itr++;
	}
}
void BioApplication::ListPresentAssociations(uint maxGeneCount) {
	string filename = reportPrefix + "-associations.txt";
	ofstream file(filename.c_str());
	map<uint, KbMetaGroup*>::iterator itr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator end = metagroups.end();

	while (itr != end) {
		if (itr->second->GetGroupCount() > 0 && !itr->second->IsDiseaseDependent()) 
			itr->second->ListAssociations(file, maxGeneCount);

		itr++;
	}
	reportLog<<setw(45)<<right<<"Association Report: "<<filename<<"\n";
}

void BioApplication::ListPopulationIDs() {
	rowset<row> rs = (sociDB.prepare << "SELECT population_label, pop_ld_comment FROM populations");
	cerr<<"Label\tComment\n";
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		cerr<<row.get<string>(0)<<"\t"<<row.get<string>(1)<<"\n";
	}
}
void BioApplication::ListGroupIDs(vector<string>& searchCriteria) {
	string pattern = "";
	string keywords = "";
	string GroupCat = " AND ";
	if (searchCriteria.size() > 0)  {
		pattern = "WHERE group_desc IN (";
		string groupFilter = "";
		string groupNameFilter = "";
		int count = searchCriteria.size();
		for (int i=0; i<count; i++ ) {
			if (groupFilter.length() > 0) {
				groupFilter+=string(" ") + GroupCat + string(" ");
				groupNameFilter+=string(" ") + GroupCat + string(" ");
			}
			groupFilter += " group_desc LIKE '%" + searchCriteria[i] + "%' ";
			groupNameFilter += " group_name LIKE '%" + searchCriteria[i] + "%' ";
		}
		pattern  = "WHERE (" + groupFilter + ") OR (" + groupNameFilter + ")";
	}
	string sql = "SELECT group_type, group_id, group_name, group_desc FROM group_type NATURAL JOIN (SELECT * FROM groups "+ pattern + ")";

	cout<<"Meta Group\tGroup ID\tName\tDescription\n";
	rowset<row> rs = (sociDB.prepare << sql);
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		cout<<row.get<string>(0)<<"\t"<<row.get<int>(1)<<"\t"<<row.get<string>(2)<<"\t"<<row.get<string>(3)<<"\n";
	}
}


void BioApplication::LoadDiseaseDependent(uint &groupID, uint &groupTypeID, ostream& os, const char *filename, const char *pop) {
	ifstream file(filename);
	KbGroup *curGroup = NULL;
	KbMetaGroup *meta = NULL;
	string regionAliases = "";
	int totalGroupCount = 0;
	bool readHeader = true;
	set<uint> groupIDsL;

	map<std::string, uint> idLookup;
	vector<string> aliases;
	GetAliasMap(idLookup);

	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Disease dependent data was not read in.\n";
		exit(1);
	}
	while (file.good()) {
		char line[4096];
		file.getline(line, 4096);
		if (readHeader) {
			readHeader = false;
			stringstream ss(line);
			string mName;
			ss>>mName;

			struct stat fStat;
			stat(filename, &fStat);
			meta = new KbMetaGroup(++groupTypeID, true, mName.c_str(), asctime(localtime(&fStat.st_mtime)) );
			cout<<"Disease Dependent ("<<filename<<" "<<mName<<")\n";
			//wrapper = new Group(groupTypeID, groupTypeID, "Wrapper", "");
			//meta->AddChild(wrapper);
		}
		else if (line[0] == '#') { }
		else if (strncmp(line, "GROUP", 5) == 0) {
			if (curGroup) {
				map<uint, KbRegion*> localRegions;
				LoadRegions(regionAliases, localRegions, os, pop);
				map<uint, KbRegion*>::iterator itr = localRegions.begin();
				map<uint, KbRegion*>::iterator end = localRegions.end();

				while (itr != end ) {
					KbRegion *region = itr->second;
					region->InsertDD(groupID);
					curGroup->Add(region);

					if (GetRegion(region->DbID()) == NULL) {
						regions[itr->first] = region;
					}
					itr++;
				}

				vector<string>::iterator aitr = aliases.begin();
				vector<string>::iterator aend = aliases.end();

				os<<"Unmatched aliases:\n";
				while (aitr != aend) {
					if (idLookup.find(*aitr) == idLookup.end()) {
						os<<"--" + *aitr<<"\n";
					}
					aitr++;
				}
				aliases.clear();

				//cout<<setw(35)<<curGroup->CommonName()<<setw(10)<<curGroup->DbID()<<setw(15)<<1<<setw(15)<<localRegions.size()<<"\n";

				regionAliases = "";
				meta->AddGroup(curGroup);
			}
			string gName, gDesc;
			stringstream ss(line +6);
			ss>>gName>>gDesc;
			//groupIDs.insert(++groupID);
			curGroup = new KbGroup(++groupID, groupID, true, gName.c_str(), gDesc.c_str());
			os<<"\n\nAdding new group ("<<groupID<<"), "<<gName<<" - "<<gDesc<<"\n";
			os<<setw(15)<<" "<<setw(20)<<"Ensembl"
				<<setw(10)<<" "<<setw(6)<<" "
				<<setw(12)<<"Start"<<setw(12)<<"Stop"
				<<setw(12)<<"Eff."<<setw(12)<<"Eff."
				<<setw(8)<<"Snp"<<"\n";
			os<<setw(15)<<"Alias"<<setw(20)<<"ID"
				<<setw(10)<<"ID"<<setw(6)<<"Chrom"
				<<setw(12)<<"Pos."<<setw(12)<<"Pos."
				<<setw(12)<<"Start"<<setw(12)<<"End"
				<<setw(8)<<"Count"<<"\n";
			os<<"-----------------------------------------------------------------------------------------------------------\n";
			totalGroupCount++;
		} else {
			if (strlen(line) > 0) {
				aliases.push_back(string(line));
				if (regionAliases.length() > 0) 
					regionAliases += ", ";
				regionAliases+="'"+string(line)+"'";
			}
		}
	}
	if (curGroup == NULL)  {
		cerr<<"Unable to load disease dependent information properly from the file, "<<filename<<". Please see the manual for instructions on the file's format.\n";
		exit(1);
	}
	if (regionAliases.length() > 0) {
		map<uint, KbRegion*> localRegions;
		//cerr<<"The aliases are: "<<regionAliases<<"\n";
		LoadRegions(regionAliases, localRegions, os, pop);
		map<uint, KbRegion*>::iterator itr = localRegions.begin();
		map<uint, KbRegion*>::iterator end = localRegions.end();
		
		while (itr != end ) {
			KbRegion *region = itr->second;

			region->InsertDD(groupID);
			curGroup->Add(region);
			regions[itr->first] = region;
			if (GetRegion(region->DbID()) == NULL) {
				regions[itr->first] = region;
			}	
			itr++;
		}
		//cout<<setw(35)<<curGroup->CommonName()<<setw(10)<<curGroup->DbID()<<setw(15)<<1<<setw(15)<<localRegions.size()<<"\n";
		meta->AddGroup(curGroup);


		vector<string>::iterator aitr = aliases.begin();
		vector<string>::iterator aend = aliases.end();
		os<<"Unmatched aliases:\n";
		while (aitr != aend) {
			if (idLookup.find(*aitr) == idLookup.end()) {
				os<<"--" + *aitr<<"\n";
			}
			aitr++;
		}
		aliases.clear();

	}
	if (totalGroupCount)
		metagroups[groupTypeID] = meta;
}

void BioApplication::AddUserDefinedGroup(const char *filename) {
	diseaseDependentFiles.push_back(filename);
}

void BioApplication::ListMetaGroups(ostream& os) {
	map<uint, KbMetaGroup*>::iterator itr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator end = metagroups.end();

	while (itr != end) {
		os<<itr->second->DbID()<<"\t"<<itr->second->Name()<<"\t"<<itr->second->GetGroupCount()<<"\n";
		itr++;
	}
}


void BioApplication::LoadGroupData(int maxSizeForActive, vector<uint>& includedGroups, ostream& os, const char *pop, const char *prefRegionNames) {
	//We need to define a few decent numbers that we can use to identify these new items with
	try {
		sociDB << "SELECT max(gene_id) FROM regions", into(maxRegionID);
		sociDB << "SELECT max(group_id) FROM groups", into(maxGroupID);
		sociDB << "SELECT max(group_type_id) FROM group_type", into(maxGroupTypeID);
	
		string groupList = "";
		vector<uint>::iterator itr = includedGroups.begin();
		vector<uint>::iterator end = includedGroups.end();
		while (itr != end) {
			if (groupList.length() > 0)
				groupList+=", ";
			groupList += Utility::ToString(*itr);
			itr++;
		}

		bool grabAll = includedGroups.size() == 0;

		uint popID = GetPopID(pop);
		rowset<row> rs = (sociDB.prepare << "SELECT group_type_id, group_type, download_date FROM group_type");
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			size_t groupID = row.get<int>(0);
			//I need to figure out how to process the date type. 
			//row.get<date> is probably the command, and I'm assuming it returns a unix timestamp...not sure though
			KbMetaGroup *meta = new KbMetaGroup(groupID, false, row.get<string>(1).c_str(), "");
			if (metagroups.find(groupID) != metagroups.end()) {
				cerr<<"Duplicate Meta Group IDd\n";
				exit(1);
			}
			metagroups[groupID] = meta;
			if (grabAll) {
				includedGroups.push_back(groupID);
			}
		}
		map<uint, string> aliasLookup;
		aliasLookup = LoadRegionAlias(prefRegionNames);
		
		cout<<"\n"<<setw(35)<<"Group"<<setw(10)<<"Group ID"<<setw(15)<<"Group-Count"<<"  "<<setw(15)<<"Gene-Count"<<"\n";
		cout<<"-----------------------------------------------------------------------------\n";
		map<uint, KbMetaGroup*>::iterator mitr = metagroups.begin();
		map<uint, KbMetaGroup*>::iterator mend = metagroups.end();
		while (mitr != mend) {
			//mitr->second->LoadGroups(sociDB, maxSizeForActive, includedGroups);
			mitr->second->LoadGroups(sociDB, maxSizeForActive, groupList.c_str(), groups);
			mitr->second->AssociateGenes(sociDB, regions, *this, aliasLookup, popID);
			//mitr->second->Prune();
			mitr++;
		}
	
		vector<string>::iterator sitr = diseaseDependentFiles.begin();
		vector<string>::iterator send = diseaseDependentFiles.end();
	
		while (sitr != send) {
			LoadDiseaseDependent(maxGroupID, maxGroupTypeID, os, (*sitr++).c_str(), pop);
		}
		InitGeneLookup(regions);

	} catch (soci_error const &e) {
		cerr<<"Unable to Load group data from database, "<<filename<<". Error: "<<e.what()<<"\n";
		exit(1);
	}
}

void BioApplication::InitGeneLookup(std::map<uint, Knowledge::KbRegion*>& regions) {
	std::map<uint, Knowledge::KbRegion*>::iterator itr = regions.begin();
	std::map<uint, Knowledge::KbRegion*>::iterator end = regions.end();

	while (itr != end) {
		GeneGeneModel::geneLookup[itr->first] = itr->second;
		itr++;
	}
}
void BioApplication::StripOptimization() {
	cout<<"Stripping optimizations from the local data-source. This is done\n"
		"to speed up certain activities (such as LD import). After that is \n"
		"completed, users should run the optimization once again.\n";
	//Test the DB connection
	sociDB << "DROP INDEX IF EXISTS group_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS group_relationships_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS group_associations_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS region_alias_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS regions_alias_aliasidx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS region_bounds_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS rs_merged_idx";
	cerr<<"*";
	sociDB << "DROP INDEX IF EXISTS regions_idx";
	cerr<<"*\n";

}

void BioApplication::PerformOptimization() {
	cout<<"Optimizing the local data-source. This operation is only necessary\n"
		"one time per data-source and could take several minutes to complete.\n";
	//Test the DB connection
	sociDB << "CREATE INDEX IF NOT EXISTS group_idx ON groups (group_type_id, group_id)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS group_relationships_idx ON group_relationships(child_id, parent_id)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS group_associations_idx ON group_associations (group_id, gene_id)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS region_alias_idx ON region_alias (region_alias_type_id, gene_id)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS region_alias_type_alias_idx ON region_alias (region_alias_type_id, alias)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS regions_alias_aliasidx ON region_alias(alias)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS region_bounds_idx ON region_bounds(gene_id, population_id)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS regions_idx ON regions (gene_id, chrom)";
	cerr<<"*";
	sociDB << "CREATE INDEX IF NOT EXISTS rs_merged_idx ON rs_merged (merged_rs_id, current_rs_id)";
	cerr<<"*\n";
}

map<uint, string> BioApplication::LoadRegionAlias(const char *filename) {
	map<uint, string> lookup;
	if (string(filename).length() == 0)
		return lookup;
	if (filename) {

		ostream *os = &cout;
		ofstream file;
		if (htmlReports) {
			string reportFilename = reportPrefix + string("-PreferredAliases.html");
			reportLog<<setw(45)<<right<<"Preferred Aliases Report: "<<reportFilename<<"\n";
			file.open(reportFilename.c_str());
			os = &file;
		}
		string aliases = "WHERE alias in ('" + Utility::FileToString(filename, "','") + "')";
		if (htmlReports) {
			(*os)<<"<HTML><HEAD><TITLE>Gene Aliases</TITLE</HEAD>\n<BODY><TABLE>\n";
			(*os)<<"<TR><TH>Alias</TH><TH>Source</TH><TH>Primary Name</TH></TR>\n";
		}
		else {
			(*os)<<"Alias (source)\tDescription\tEntrez Gene Reference\n";
		}
		rowset<row> rs = (sociDB.prepare << ("SELECT gene_id, alias, region_alias_type_desc, primary_name FROM region_alias_type NATURAL JOIN region_alias NATURAL JOIN regions " + aliases));
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			size_t geneID = row.get<int>(0);
			string alias = row.get<string>(1);
			//string desc = row.get<string>(2);
			string sourceDesc = row.get<string>(2);
			string name = row.get<string>(4);
			if (htmlReports)
				(*os)<<"<TR><TD><A HREF='http://www.ncbi.nlm.nih.gov/gene/"<<geneID<<"'>"<<alias<<"</A></TD><TD>"
					<<sourceDesc<<"</TD><TD>"<<name<<"</TD></TR>\n";
			else
				(*os)<<alias<<" ("<<sourceDesc<<")\t"<<name<<"\thttp://www.ncbi.nlm.nih.gov/gene/"<<geneID<<"\n";
			lookup[geneID] = alias;
		}
		if (htmlReports)
			(*os)<<"</TABLE></BODY></HTML>\n";
		//cerr<<"ASDFASDFASDFASDFASDF "<<lookup.size()<<"\n";
	}

	std::map<uint, Knowledge::KbRegion*>::iterator rItr = regions.begin();
	std::map<uint, Knowledge::KbRegion*>::iterator rEnd = regions.end();
	while (rItr != rEnd) {
		uint regionID										= rItr->first;
		if (lookup.find(regionID) == lookup.end())
			lookup[regionID]								= rItr->second->CommonName();
		rItr++;
	}

	return lookup;
}

void BioApplication::ListGenes(std::ostream& os) {
	//SetRegionBoundaries("NO-LD");

	std::map<uint, Knowledge::KbRegion*>::iterator itr = regions.begin();
	std::map<uint, Knowledge::KbRegion*>::iterator end = regions.end();
	os<<"Name\tChr\tStart\tEnd\tStrand\tEnsembl ID\tEff. Start\tEff. Endn";
	while (itr != end) {
		if (itr->second->SnpCount() > 0) {
			uint start, stop;
			itr->second->GetOrigBounds(start, stop);
			uint effStart, effStop;
			itr->second->GetBounds(effStart, effStop);
			//For now, we don't have strand information, so we will just pass an empty cell
			os<<itr->second->CommonName()<<"\t"<<itr->second->Chromosome()<<"\t"<<start<<"\t"<<stop<<"\t"<<"\t"<<itr->second->Name()<<"\t"<<effStart<<"\t"<<effStop<<"\n";
		}
		itr++;
	}

}
void BioApplication::BuildSnpToGene(SnpToGeneManager& snpsToGenes, std::vector<std::string>& geneList) {
	std::map<uint, Knowledge::KbRegion*>::iterator regionItr = regions.begin();
	std::map<uint, Knowledge::KbRegion*>::iterator regionEnd = regions.end();

	while (regionItr != regionEnd) {
		std::set<uint> snps = regionItr->second->SNPs();
		std::set<uint>::iterator itr = snps.begin();
		std::set<uint>::iterator end = snps.end();
		while (itr != end) {
			uint rsID = *itr++;
			//cerr<<regionItr->second->Chromosome()<<","<<regionItr->second->CommonName()<<",rs"<<rsID<<"\n";
			snpsToGenes.AddGene(regionItr->second->Chromosome().c_str(), rsID, regionItr->second);
			//os<<regionItr->second->Chromosome()<<","<<regionItr->second->CommonName()<<",rs"<<*itr++<<"\n";
		}
		regionItr++;
	}
}
void BioApplication::ReportSNPsInGenes(std::vector<std::string>& genelist, ostream& os) {
	std::map<uint, Knowledge::KbRegion*>::iterator regionItr = regions.begin();
	std::map<uint, Knowledge::KbRegion*>::iterator regionEnd = regions.end();

	while (regionItr != regionEnd) {
		std::set<uint> snps = regionItr->second->SNPs();
		std::set<uint>::iterator itr = snps.begin();
		std::set<uint>::iterator end = snps.end();
		while (itr != end) {
			os<<regionItr->second->Chromosome()<<","<<regionItr->second->CommonName()<<",rs"<<*itr++<<"\n";
		}
		regionItr++;
	}
}

void BioApplication::DetailCoverage(std::vector<std::string>& genelist, std::vector<std::string>& snpFiles, bool detailedCoverage) {
	string filename = reportPrefix + ".gene-coverage";
	if (htmlReports) 
		filename = reportPrefix + "-gene-coverage.html";

	ofstream os(filename.c_str());
	vector< set<uint> >snpSets;

	if (htmlReports) {
		os<<"<HTML>\n<HEAD>\n<TITLE>Gene Coverage</TITLE>\n";
		os<<"<script language=\"javascript\">\n<!--\n	var state = 'none';\n	function showhide(layer_ref) {\n\n	if (state == 'block') {\n		state = 'none';\n	} else {\n		state = 'block';\n	}\n	if (document.all) { //IS IE 4 or 5 (or 6 beta)\n		eval( \"document.all.\" + layer_ref + \".style.display = state\");\n	}\n	if (document.layers) { //IS NETSCAPE 4 or below\n		document.layers[layer_ref].display = state;\n	}\n	if (document.getElementById &&!document.all) {\n		hza = document.getElementById(layer_ref);\n		hza.style.display = state;\n	}\n}\n//-->\n</script> \n</HEAD>\n<BODY>\n<TABLE CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES>\n";
		if (detailedCoverage)
			os<<"<TR bgcolor='#F3EFE0'><TH>Gene</TH><TH>Ensembl ID</TH><TH>Chromosome</TH><TH>Begin(kB)</TH><TH>End(kB)<TH>Effective Begin</TH><TH>Effective End</TH></TH><TH>Total</TH>";
		else
			os<<"<TR bgcolor='#F3EFE0'><TH>Gene</TH><TH>Ensembl ID</TH><TH>Total</TH>";
	}
	else
		if (detailedCoverage)
			os<<"Gene\tEnsembl_id\tChromosome\tBegin(kB)\tEnd(kB)\tEffective Begin(kB)\tEffective End(kB)\tTotal\t";
		else
			os<<"Gene\tEnsembl_id\tTotal\t";

	for (size_t i=0; i<snpFiles.size(); i++) {
		set<uint> snpSet;
		ifstream snpFile(snpFiles[i].c_str());
		while (snpFile.good()){ 
			size_t snp = 0;
			snpFile>>snp;
			if (snp > 0)
				GetSNPs(snp, snpSet);
		}
		if (htmlReports) {
			os<<"<TH>"<<snpFiles[i]<<"("<<snpSet.size()<<")</TH>";
			if (detailedCoverage)
				os<<"<TH>SNPs</TH>";
		} else {
			os<<snpFiles[i]<<"("<<snpSet.size()<<")\t";
			if (detailedCoverage)
				os<<"SNPs\t";
		}
		snpSets.push_back(snpSet);
	}
	if (htmlReports)
		os<<"</TR>\n";
	else
		os<<"\n";
	map<string, KbRegion*> genes;
	vector<string>::iterator itr = genelist.begin();
	vector<string>::iterator end = genelist.end();
	while (itr != end) {
		size_t geneID = 0;
		sociDB << "SELECT gene_id FROM region_alias WHERE alias=:gene", use(*itr, "gene"), into(geneID);
		if (htmlReports)
			os<<"<TR><TD>"<<*itr<<"</TD>";
		else
			os<<*itr<<"\t";
		if (geneID > 0 && regions.find(geneID) != regions.end()) {
			KbRegion *region = regions[geneID];
			if (htmlReports)
				os<<"<TD><A HREF='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g="<<region->Name()<<"'>"<<region->Name()<<"</A></TD>";
			else
				os<<region->Name();
			if (htmlReports) {
				if (detailedCoverage)
					os<<"<TD>"<<region->Chromosome()<<"</TD><TD>"<<region->OriginalStart()<<"</TD><TD>"<<region->OriginalEnd()<<"</TD><TD>"<<region->Start()<<"</TD><TD>"<<region->End()<<"</TD>";
				os<<"<TD>"<<region->SnpCount()<<"</TD>";
			} else {
				if (detailedCoverage)
					os<<"\t"<<region->Chromosome()<<"\t"<<region->OriginalStart()<<"\t"<<region->OriginalEnd()<<"\t"<<region->Start()<<"\t"<<region->End();
				os<<"\t"<<region->SnpCount()<<"\t";
			}
			for (size_t i=0; i<snpSets.size(); i++) {
				set<SNP_Details> snps; 
				region->GetSnpCoverage(snpSets[i], snps);
		
				if (htmlReports)
					os<<"<TD>"<<snps.size()<<"</TD>";
				else
					os<<snps.size()<<"\t";
				if (detailedCoverage) {
					set<SNP_Details>::iterator itr = snps.begin();
					set<SNP_Details>::iterator end = snps.end();
					bool isFirst=true;
					if (htmlReports) 
						os<<"<TD>";
					int count = 0;
					while (itr != end) {
						if (htmlReports) {
							if (count++ == 10)
								os<<"<DIV id=\""<<region->Name()<<"."<<i<<"\" style=\"display: none;\">";
							os<<"\n\t\t<A HREF='http://www.ensembl.org/Homo_sapiens/Variation/Summary?source=dbSNP;v=rs"<<itr->rsID<<"'>rs"<<itr->rsID<<"</A> ";
						}
						else {
							if (!isFirst) {
								os<<" ";
							}
							os<<"rs"<<itr->rsID<<" ( "<<itr->position<<" )";
						}
						itr++;
						isFirst = false;
					}
					if (htmlReports) {
						if (count>9) 
							os<<"</DIV><A HREF=\"#"<<region->Name()<<"\" onclick=\"showhide('"<<region->Name()<<"."<<i<<"');\">...</A>";
						os<<"</TD>";
					}
					else
						os<<"\t";
				}
			}
		}
		if (htmlReports)
			os<<"</TR>";
		os<<"\n";
		itr++;
	}
	if (htmlReports)
		os<<"</TABLE>\n<P>*Boundaries for genes are based on: "<<populationDesc<<".\n</BODY></HTML>";
	else
		os<<"\n*Boundaries for genes are based on: "<<populationDesc<<"\n";
	reportLog<<setw(45)<<right<<"Gene Coverage Report: "<<filename<<"\n";
}

void BioApplication::SummarizeModelCounts(int maxGeneCount) {
	GeneGeneModelArchive geneModels;

	map<uint, KbMetaGroup*>::iterator itr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator end = metagroups.end();
	uint modelCount = 0;

	while (itr != end) {
		uint localCount = itr->second->GenerateGeneGeneModels(geneModels, maxGeneCount, cout);
		modelCount += localCount;
		cout<<setw(40)<<itr->second->Name()<<" Gene Model Count: "<<localCount<<"\n";
		itr++;
	}
	cout<<"\nTotal Gene-Gene Model Count: "<<modelCount<<"\n";
	map<uint, uint> modelCounts;

	std::map<uint, Region*> *r = (std::map<uint, Region*>*)&regions;
	geneModels.SummarizeModelCounts(modelCounts, *r);


	map<uint, uint>::iterator impl = modelCounts.begin();
	map<uint, uint>::iterator implend = modelCounts.end();
	while (impl != implend) {
		modelCount += (uint) (impl++->second);
	}
	impl = modelCounts.begin();
	
	cerr<<"Counts by group pairings:\n";
	while (impl != implend) {
		cerr<<impl->first<<"\t"<<impl->second<<"\t"<<((float)impl->second/(float)modelCount*100.0)<<"%\n";
		impl++;
	}
	cerr<<"Total Model Count: "<<modelCount<<"\n";
}

/*
struct Write {
	ostream &os;
	map<uint, KbMetaGroup*>& metaGroups;
	map<uint, KbRegion*>& regions;

	Write(ostream& os, map<uint, KbMetaGroup*>& metaGroups, map<uint, KbRegion*>& regions) : os(os), metaGroups(metaGroups), regions(regions) { }

	void operator()(const Model& model) {
		model.Write(os, metaGroups, regions);
	}
};
*/

void BioApplication::SetVariationFilename(const char *filename) {
	sociDB<<"UPDATE versions SET version=:var WHERE element='variations'", use(string(filename));
}


void BioApplication::SnpReport(ostream& os, ostream& failedSnps, vector<uint>& snps, bool writeHTML) {
	Reporting::SnpToGeneMapping fn(os, failedSnps, writeHTML);
	fn(regions, snps, *this);
}

void BioApplication::ProduceModels(GeneGeneModelArchive& geneModels, ostream& os, int maxGeneCount) {
	map<uint, KbMetaGroup*>::iterator itr = metagroups.begin();
	map<uint, KbMetaGroup*>::iterator end = metagroups.end();
	uint modelCount = 0;
	os<<"----------------------------------------Gene-Gene Models------------------------------------------------------------\n";
	os<<setw(35)<<"Gene"<<setw(8)<<"SNP"<<setw(35)<<"Gene"<<setw(8)<<"SNP"<<setw(10)<<"Impl."<<setw(10)<<"Models"<<"\tGroups\n";
	os<<setw(35)<<"Name"<<setw(8)<<"Count"<<setw(35)<<"Name"<<setw(8)<<"Count"<<setw(10)<<"Index"<<setw(10)<<"Count"<<"\tDI,DD\n";
	os<<"--------------------------------------------------------------------------------------------------------------------\n";

	while (itr != end) {
		modelCount += itr->second->GenerateGeneGeneModels(geneModels, maxGeneCount, os);
		itr++;
	}
	cout<<"\nTotal Gene-Gene Model Count: "<<modelCount<<"\n";
	
	//geneModels.GenerateModels(repo, regions, os, maxGeneCount, 0);
}

int BioApplication::GetPopID(const char *pop) {
	int popID = -1;
	sociDB<<"SELECT population_id FROM populations WHERE population_label = :pop", use(string(pop)), into(popID);
	return popID;
}

string BioApplication::GetPopulationDesc(const char *pop) {
	string desc;
	sociDB<<"SELECT pop_ld_comment FROM populations WHERE population_label = :pop", use(string(pop)), into(desc);
	return desc;
}
uint BioApplication::LoadRegions(const string& geneList, const char *pop) {
	stringstream ss;
	populationDesc = GetPopulationDesc(pop);
	return LoadRegions(geneList, regions, ss, pop);
}

void BioApplication::LoadRegionAliases(const string& geneList) {
	string filename = string(reportPrefix) + ".aliases";

	if (htmlReports) 
		filename = reportPrefix + string("-aliases.html");

	ofstream file(filename.c_str());
	if (htmlReports) {
		file<<"<HTML><HEAD><TITLE>Region Aliases</TITLE></HEAD>\n<BODY><TABLE>\n";	
		file<<"<TR><TH>Gene Alias</TH><TH>Source</TH><TH>Ensembl ID</TH><TH>Description</TH></TR>\n";
	}
	else
		file<<"\nGene Aliases: \nGene Alias\tSource\tEnsembl\tDescription\n";
	map<uint, KbRegion*>::iterator notFound = regions.end();

	string geneFilter = "";
	if (geneList != "ALL" and geneList != "''")
		geneFilter					= " AND alias IN (" + geneList + ") ";

	//cerr<<"SELECT gene_id, alias, alias_desc, region_alias_type_desc, ensembl_id FROM region_alias_type NATURAL JOIN region_alias NATURAL JOIN regions WHERE region_alias_type_id=1300 " + geneFilter + " ORDER BY ensembl_id, alias DESC\n";
	rowset<row> rs = (sociDB.prepare << "SELECT gene_id, alias, description, region_alias_type_desc, primary_name FROM region_alias_type NATURAL JOIN region_alias NATURAL JOIN regions WHERE region_alias_type_id IN (1300, 2000, 2200) " + geneFilter + " ORDER BY gene_id, alias DESC");
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		size_t geneID = row.get<int>(0);
		string alias = row.get<string>(1);	
		string desc = row.get<string>(2);
		string sourceDesc = row.get<string>(3);
		string name = row.get<string>(4);
		if (regions.find(geneID) != notFound) {
			regions[geneID]->AddAlias(alias.c_str());

			if (htmlReports)
				file<<"<TR><TD>"<<alias<<"</TD><TD>"<<sourceDesc<<"</TD><TD><A HREF='http://www.ncbi.nlm.nih.gov/gene/"<<geneID<<"'>"<<name<<"</A></TD><TD>"<<desc<<"</TD></TR>\n";
			else
				file<<alias<<" ("<<sourceDesc<<")\t"<<desc<<"\thttp://www.ncbi.nlm.nih.gov/gene/"<<geneID<<"\n";
		}
	}
	cerr<<regions.size()<<" regions loaded.\n";
	reportLog<<setw(45)<<right<<"Alias Report: "<<filename<<"\n";
}

void BioApplication::SetRegionBoundaries(const char *pop, map<uint, KbRegion*>& regions) {
	int popID = GetPopID(pop);
	if (popID > 0)
		geneExtension		= 0;
	string sql = "SELECT gene_id, chrom, start, end  FROM regions NATURAL JOIN region_bounds WHERE population_id="+Utility::ToString(popID);

	//cerr<<sql<<"\n";
	//uint geneCount = 0;

	try {
		rowset<row> rs = (sociDB.prepare << sql);
		//uint snpCount = 0;

		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			uint gene_id = row.get<int>(0);
			uint start = row.get<int>(2);
			if (start < geneExtension)
				start			= 1;
			else
				start			= start - geneExtension;

			uint stop		= row.get<int>(3) + geneExtension;

			if (regions.find(gene_id) != regions.end())
				regions[gene_id]->SetBounds(start, stop);
		}
	}
	catch (exception const &e) {
		cerr<<sql<<"\n";
		cerr<<"Unable to read region data from the database. DB Error: "<<e.what()<<"\n";
	}

}
uint BioApplication::LoadVariations(set<uint>& snps, const char *mergedReport, bool doReport) {
	string filename;
	string path, basename, extension;
	Utility::SplitIntoComponents(dbFilename.c_str(), path, basename, extension);
	sociDB<<"SELECT version FROM versions WHERE element = 'variations'", into(filename);

	if (filename.length() == 0) {
		cerr<<"There is a problem with the database. It doesn't have a link to a variations file.";
		abort();
	}

	if (!Utility::FileExists(filename.c_str())) {
		filename = path + "/" + filename;
	}
	return InitSNPs(snps, filename.c_str(), mergedReport, doReport);
}
void BioApplication::GetAliasMap(map<string, uint>& idLookup) {
	string sql						= "SELECT gene_id, alias FROM region_alias WHERE region_alias_type_id IN (1300,2000,2200)";
	try {
		rowset<row> rs = (sociDB.prepare << sql);
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row			= *itr;
			idLookup[row.get<string>(1)] = row.get<int>(0);
		}
	} catch (exception const &e) {
		cerr<<sql<<"\n";
		cerr<<"An error was encountered during production of the alias lookup table. DB Error: "<<e.what()<<"\n";
	}
}

uint BioApplication::LoadRegions(const string& geneList, map<uint, KbRegion*>& regions, ostream& os, const char *pop) {
	int popID = GetPopID(pop);
	if (popID > 0)
		geneExtension				= 0;
	string geneFilter = " WHERE ";
	if (geneList != "ALL" && geneList != "''")
		geneFilter					= "JOIN ( SELECT * FROM region_alias WHERE alias IN (" + geneList + ") ) a ON (a.gene_id=regions.gene_id) AND ";
	//string sql = "SELECT m.gene_id, ensembl_id, chrom, start, end, description, a.alias FROM (SELECT * FROM regions NATURAL JOIN region_bounds WHERE population_id="+Utility::ToString(popID)+") m INNER JOIN  (SELECT * FROM region_alias WHERE region_alias_type_id IN (1300) " + geneFilter + ") a ON (m.gene_id=a.gene_id) ORDER BY alias DESC";
	string sql						= "SELECT regions.gene_id, primary_name, chrom, region_bounds.start, region_bounds.end, description, orig.start AS orig_start, orig.end AS orig_end  FROM regions NATURAL JOIN region_bounds JOIN (SELECT * FROM region_bounds WHERE population_id=0) AS orig ON (regions.gene_id=orig.gene_id) " + geneFilter +" region_bounds.population_id=" + Utility::ToString(popID);

	uint geneCount = 0;

	try {
		rowset<row> rs = (sociDB.prepare << sql);
		uint snpCount = 0;
	
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row			= *itr;
			uint gene_id			= row.get<int>(0);
			string name			= row.get<string>(1);
			geneCount++;
			KbRegion *region		= NULL;

			region = GetRegion(gene_id);
			if (region == NULL) {
				if (regions.find(gene_id) == regions.end()) {
					uint start		= row.get<int>(3);
					if (start < geneExtension)
						start			= 1;
					else
						start			= start - geneExtension;

					uint stop		= row.get<int>(4) + geneExtension;
					string chrom	= row.get<string>(2);
					string desc		= "";
					uint origStart = row.get<int>(6);
					uint origEnd   = row.get<int>(7);
					if (row.get_indicator(5) != i_null)
						desc = row.get<string>(5);
					region = new KbRegion(gene_id, start, stop, chrom.c_str(), name.c_str(),desc.c_str(), this);
					region->SetOrigBounds(origStart, origEnd);
					snpCount = region->AssociateSNPs();
					regions[gene_id] = region;
				}
				else {
					region			= regions[gene_id];
				}
			}
			regions[gene_id]		= region;
			//uint start, stop;
			//region->GetBounds(start, stop);
		}
	} catch (exception const &e) {
		cerr<<sql<<"\n";
		cerr<<"Unable to read region data from the database. DB Error: "<<e.what()<<"\n";
	}
	SetRegionBoundaries(pop, regions);
	LoadRegionAliases(geneList);
	map<uint, KbRegion*>::iterator itr = regions.begin();
	map<uint, KbRegion*>::iterator end = regions.end();

	while (itr != end) {
		KbRegion *r = itr++->second;
		os<<setw(15)<<r->CommonName()<<setw(20)<<r->Name()<<setw(10)<<r->DbID()<<setw(6)<<r->Chromosome()
			<<setw(12)<<r->Start()<<setw(12)<<r->End()<<setw(12)<<r->OriginalStart()<<setw(12)<<r->OriginalEnd()<<setw(5)<<r->SnpCount()<<"\n";
	}
	return geneCount;
}

void BioApplication::LoadConversion() {
	rsConversion.clear();
	string sql = "SELECT merged_rs_id, current_rs_id FROM rs_merged";

	try {
		rowset<row> rs = (sociDB.prepare<<sql);
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			uint rsid						= row.get<int>(0);
			uint newrs						= row.get<int>(1);

			rsConversion[rsid]			= newrs;
		}
		ofstream file("conversions.csv");
		map<uint, uint>::iterator itr = rsConversion.begin();
		map<uint, uint>::iterator end = rsConversion.end();

		while (itr != end) {
			file<<itr->first<<"\t"<<itr->second<<"\n";
			itr++;
		}
	} catch (exception const &e) {
		cerr<<"Unable to load merge information from dbSNP. DB Error: "<<e.what()<<"\n"<<sql<<"\n";
	}


}


void BioApplication::InitBiofilter(const char *dbFilename, bool reportVersion) {
	//Test the DB connection
	this->dbFilename							= dbFilename;
	if (!Utility::FileExists(dbFilename)) {
		cerr<<"The database, "<<dbFilename<<", could not be found. Unable to continue.\n";
		exit(1);
	}
	try {
		string cnxParam = "dbname="+string(dbFilename)+" timeout=2500";
		sociDB.open(soci::sqlite3, cnxParam.c_str());
		string dbSnp, ensembl, hapmap, build, variations;
		sociDB<<"SELECT version FROM versions WHERE element='ncbi'", into(dbSnp);
		sociDB<<"SELECT version FROM versions WHERE element='ensembl'", into(ensembl);
		sociDB<<"SELECT version FROM versions WHERE element='hapmap'", into(hapmap);
		sociDB<<"SELECT version FROM versions WHERE element='variations'", into(variations);
		sociDB<<"SELECT version FROM versions WHERE element='build'", into(build);
		this->varVersion				= atoi(variations.c_str());

		if (reportVersion) {
			cerr<<"\n------------------------- Dependency Versions ----------\n";
			cerr<<setw(35)<<right<<"dbSNP: "<<dbSnp<<"\n";
			cerr<<setw(35)<<right<<"Ensembl: "<<ensembl<<"\n";
			cerr<<setw(35)<<right<<"Hap Map LD: "<<hapmap<<"\n";
			cerr<<setw(35)<<right<<"Variation Filename: "<<variations<<"\n";
			cerr<<setw(35)<<right<<"Genome Build: "<<build<<"\n";

		}
		
	} catch (soci::soci_error const &e) {
		cerr<<"Problems were encountered trying to open the database, "<<dbFilename<<". Error: "<<e.what()<<"\n";
	}


	try {
		rowset<row> rs = (sociDB.prepare << "SELECT id, role FROM snp_role");

		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			uint id = row.get<int>(0);
			string name = row.get<string>(1);

			Chromosome::roleDescription[id] = name;
		}
	} catch (soci::soci_error const &e) {
		cerr<<"An error was encountered trying to load the role data. SNP role information will not be available.\n";
	}
}

}
