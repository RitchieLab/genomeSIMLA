//
// C++ Implementation: kbmetagroup.cpp
//
// Description: Primary functionality associated with a single biofilter knowledge source
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//

#include "kbmetagroup.h"
#include <iomanip>

using namespace soci;
using namespace std;
using namespace Utility;

namespace Biofilter {
namespace Knowledge {

KbMetaGroup::KbMetaGroup(uint id, bool isDiseaseDependent, const char *name, const char*date) : KbEntity(id, name), date(date), isDiseaseDependent(isDiseaseDependent) { }
KbMetaGroup::KbMetaGroup(uint id, bool isDiseaseDependent) : KbEntity(id, ""), date(""), isDiseaseDependent(isDiseaseDependent) { }
KbMetaGroup::~KbMetaGroup() {
	groups.clear();
}

void KbMetaGroup::AddGroup(KbGroup *group) {
	groups[group->DbID()] = group;
}
KbMetaGroup::iterator KbMetaGroup::begin() {
	return groups.begin();
}

KbMetaGroup::iterator KbMetaGroup::end() {
	return groups.end();
}

uint KbMetaGroup::GetGroupCount() {
	return groups.size();
}

bool KbMetaGroup::IsDiseaseDependent() {
	return isDiseaseDependent;
}
/**
 * @brief The following is an attempt to speed up the group loading process by performing one large query for each metagroup, and building associations from those in memory.
 */
uint KbMetaGroup::LoadGroups(soci::session& sociDB, int maxThesh, vector<uint>& rootSource, GroupLookup& grpLookup) {
	string sqlString = "SELECT group_id, group_name, group_desc FROM groups WHERE group_type_id=?";
	rowset<row> rs = (sociDB.prepare << "SELECT group_id, group_name, group_desc FROM groups WHERE group_type_id=?", use(dbID));
	vector<uint> roots = rootSource;
	//If it's empty, we want to grab everything
	if (roots.size() == 0)
		roots.push_back(dbID);
	for (rowset<row>::const_iterator itr=rs.begin(); itr!=rs.end(); ++itr) {
		row const& row = *itr;
		uint groupID = row.get<int>(0);
		KbGroup *group = new KbGroup(groupID, DbID(), isDiseaseDependent, row.get<string>(1).c_str(), row.get<string>(2).c_str());
		groups[groupID] = group;
		grpLookup[groupID] = group;
	}
	vector<uint>::iterator itr = roots.begin();
	vector<uint>::iterator end = roots.end();
	map<uint, KbGroup*>::iterator gEnd = groups.end();
	//All groups are loaded. Now, we need to mark those that fall within the requested groups and traverse any relationships
	//If the metagroup ID is in the set, we already have what we wanted
	if (find(itr, end, DbID()) != end) {
		while (itr != end) {
			//This is one of the requested nodes. increment the count and traverse
			if (groups.find(*itr) != gEnd) {
				KbGroup *g = groups[*itr];
				set<uint> parents;
				g->LoadAssociations(sociDB, groups, parents);
			}
			itr++;
		}
	}
	return groups.size();
}
/**
 * we are assuming that there can be no cross-relationships between meta-groups. This means
 * that reactome and other types that have a real tree relationship can not be designated
 * as having separate groups. What we can do is relate them in the tree in a way that
 * segregates the appropriately.
 */
uint KbMetaGroup::LoadGroups(soci::session& sociDB, int maxThresh, const char *inc, GroupLookup& grpLookup) {
	string inClause = "";
	if (strlen(inc) > 0)
		inClause = string("group_id IN (") + string(inc) + string(") AND");
	string sqlString = "SELECT group_id, group_name, group_desc FROM groups WHERE " + inClause + " group_type_id=?";

	rowset<row> rs = (sociDB.prepare << sqlString.c_str(), use(dbID));
	//rowset<row> rs = (sociDB.prepare << "SELECT group_id, group_name, group_desc FROM groups WHERE group_id IN '?' AND group_type_id=?", use(inc), use(id));
	uint count = 0;
	//Load each of the meta-groups
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		size_t groupID = row.get<int>(0);
		if (groups.find(groupID) != groups.end()) {
			cerr<<"!!!!!!!!!!!!!!!! Duplicate group ID: "<<groupID<<"\n";
			exit(1);
		}
		KbGroup *group = new KbGroup(groupID, DbID(), isDiseaseDependent, row.get<string>(1).c_str(), row.get<string>(2).c_str());
		//group->AddDdEntry(id);
		groups[groupID] = group;
		grpLookup[groupID] = group;
		count++;
	}

	
	//Lets load group relationships. This had to wait until the groups were loaded
	map<uint, KbGroup*>::iterator itr = groups.begin();
	map<uint, KbGroup*>::iterator end = groups.end();
	while (itr != end) {
		//We are going to expect the group to load any children that haven't
		//already been loaded, and pass an empty "parent" list. This will grow
		//as you drill further down
		set<uint> parents;
		itr->second->LoadAssociations(sociDB, groups, parents);
		itr++;
	}
	MarkAsProcessed(true);
	return count;
}


void KbMetaGroup::GraphAssociations(ostream& os, uint maxGeneCount) {
	map<uint, KbGroup *>::iterator itr = groups.begin();
	map<uint, KbGroup *>::iterator end = groups.end();


	//os<<"Associations ("<<name<<"):\n";
	//os<<"\t"<<DbID()<<"[label=\""<<CommonName()<<"\"];\n";
	while (itr != end) {
		stringstream ss;
		if (itr->second->GraphAssociations(ss, maxGeneCount) > 0) {
			//if (DbID() != itr->second->DbID())
			//	os<<"\t"<<itr->second->DbID()<<"->"<<DbID()<<"\n";
			os<<ss.str();
		}
		itr++;
	}
}
void KbMetaGroup::ListAssociations(ostream& os, uint maxGeneCount) {
	map<uint, KbGroup *>::iterator itr = groups.begin();
	map<uint, KbGroup *>::iterator end = groups.end();

	os<<"Associations ("<<name<<"):\n";
	while (itr != end) {
		stringstream ss;
		if (itr->second->ListAssociations(1, ss, maxGeneCount) > 0)
			os<<ss.str();
		itr++;
	}
}



uint KbMetaGroup::GenerateGeneGeneModels(GeneGeneModelArchive& archive, int maxGeneCount, ostream& os) {
	uint modelCount = 0;

	map<uint, KbGroup *>::iterator itr = groups.begin();
	map<uint, KbGroup *>::iterator end = groups.end();

	while (itr != end)
		modelCount+= itr++->second->GenerateGeneGeneModels(archive, DbID(),  maxGeneCount, os);

	return modelCount;
}



uint KbMetaGroup::AssociateGenes(soci::session& sociDB, map<uint,KbRegion*>& genes, SnpManager& snps, map<uint, string>& aliases, uint popID) {
	string groupIDs = "";
	map<uint, KbGroup*>::iterator itr = groups.begin();
	map<uint, KbGroup*>::iterator end = groups.end();
	map<uint, string>::iterator notFound = aliases.end();
	if (groups.size() == 0)
		return 0;

	while (itr != end) {
		if (groupIDs.length() > 0)
			groupIDs+=",";
		groupIDs+=ToString(itr->first);
		itr++;
	}
	string sql = "SELECT group_id, regions.gene_id, primary_name, chrom, start, end, description FROM (SELECT * FROM region_bounds WHERE population_id=" + ToString(popID) + ") NATURAL JOIN regions INNER JOIN group_associations ON (regions.gene_id=group_associations.gene_id) WHERE group_id IN (" + groupIDs +") ORDER BY group_id";
//cerr<<sql<<"\n";
	rowset<row> rs = (sociDB.prepare<<sql);
	//rowset<row> rs = (sociDB.prepare << "SELECT group_id, regions.gene_id, ensembl_id, chrom, start, end, description FROM (SELECT * FROM region_bounds WHERE population_id=:popID) NATURAL JOIN regions INNER JOIN group_associations ON (regions.gene_id=group_associations.gene_id) WHERE group_id in [:groups]", use(popID), use(groupIDs));
	set<uint> geneIDs;
	KbGroup *group = NULL;
	uint lastID = 0;
	cout<<setw(35)<<name<<setw(10)<<DbID()<<setw(15)<<groups.size();cout.flush();

	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint groupID = row.get<int>(0);
		uint gene_id = row.get<int>(1);
		geneIDs.insert(gene_id);
		KbRegion *region = NULL;
		if (genes.find(gene_id) == genes.end()) {
			uint start = row.get<int>(4);
			uint stop = row.get<int>(5);
			string chrom = row.get<string>(3);
			string name = row.get<string>(2);
			string desc = "";		//row.get<string>(5);
			region = new KbRegion(gene_id, start, stop, chrom.c_str(), name.c_str(),desc.c_str(), &snps);
			if (aliases.find(gene_id) != notFound)
				region->Name(aliases[gene_id].c_str());
			region->AssociateSNPs();
			genes[gene_id] = (KbRegion*)region;
		}
		if (region == NULL)
			region = genes[gene_id];
		if (lastID != groupID)
			group = (KbGroup*)groups[groupID];
		group->AssociateGene(region);
	}
	if (geneIDs.size() > 0)
		cout<<setw(15)<<geneIDs.size();
	cout<<"\n";
	return geneIDs.size();

}


	}
}
