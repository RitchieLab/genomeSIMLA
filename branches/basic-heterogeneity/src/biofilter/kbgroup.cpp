//
// C++ Implementation: kbgroup.cpp
//
// Description: Primary functionality associated with a "group" from one of the knowledge sources
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
// C++ Implementation: group
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kbgroup.h"

using namespace soci;
using namespace std;
using namespace Utility;

namespace Biofilter {

namespace Knowledge {

KbGroup::DiseaseDependentRelationshipType KbGroup::DiseaseDependentRelationship = AllModels;
bool KbGroup::CollapseAssociationReport		= false;
	

KbGroup::KbGroup(uint id, uint typeID, bool isDiseaseDependent, const char *name, const char *desc)
		: KbEntity(id, name, desc), groupType(typeID), isDiseaseDependent(isDiseaseDependent) {
	//cout<<"KbGroup: "<<id<<" "<<typeID<<" "<<isDiseaseDependent<<" "<<name<<" "<<desc<<"\n";
}
	

	
	
KbGroup::~KbGroup()	{
	children.clear();
	regions.clear();
}


void KbGroup::Add(KbRegion* region) {
	regions.insert(region);
//	region->Retain();
}

void KbGroup::Add(KbGroup* group) {
	children.insert(group);
//	group->Retain();
}


uint KbGroup::GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive, uint pairID, uint maxGeneCount, ostream& os) {
	//First, pass the structures to all children, in case they have the potential to generate models as well
	uint modelCount = 0;

	set<KbRegion*> regions;
	GetAllRegions(regions);

	if (DiseaseDependentRelationship == DD_GroupLevel) {
		set<KbRegion*>::iterator itr = regions.begin();
		set<KbRegion*>::iterator end = regions.end();
		set<uint> ddGroups;
		while (itr != end) {
			KbRegion *region = *itr++;
			set<uint> localDD = region->DdGroups();
			ddGroups.insert(localDD.begin(), localDD.end());
		}
		if (ddGroups.size() == 0)
			return 0;
	}

	if (pairID == (uint)-1 || regions.size() <= maxGeneCount) {
		set<KbRegion*>::iterator itr = regions.begin();
		set<KbRegion*>::iterator end = regions.end();

		//modelCount = regions.size();

		while (itr != end) {
			KbRegion *cur = *itr;
			set<KbRegion*>::iterator remainder = ++itr;

			while (remainder != end) {
				KbRegion *other = *remainder++;
				GeneGeneModel model(cur, other, pairID);
				if ((DiseaseDependentRelationship != DD_Only) || (DiseaseDependentRelationship == DD_Only && model.DdGroups().size() > 0)) {
					geneArchive.Insert(model);
					model.WriteSummary(os);
					modelCount++;
				}
			}
		}
	} else {
		set<KbGroup*>::iterator itr = children.begin();
		set<KbGroup*>::iterator end = children.end();

		while (itr != end)
			modelCount += (*itr++)->GenerateGeneGeneModels(geneArchive, pairID, maxGeneCount, os);
	}
	return modelCount;
}

void KbGroup::GetAllGroups(std::set<KbGroup*>& allgroups) const {
	allgroups.insert(children.begin(), children.end());
	set<KbGroup*>::iterator itr = children.begin();
	set<KbGroup*>::iterator end = children.end();
	while (itr != end)
		(*itr++)->GetAllGroups(allgroups);
}

void KbGroup::GetAllRegions(std::set<KbRegion*>& allregions) {
	set<KbGroup*>::iterator itr = children.begin();
	set<KbGroup*>::iterator end = children.end();

	set<KbRegion*>::iterator rItr = regions.begin();
	set<KbRegion*>::iterator rEnd = regions.end();
	while (rItr != rEnd) {
		if ((*rItr)->SnpCount() > 0)
			allregions.insert(*rItr);
		rItr++;
	}

	while (itr != end)
		(*itr++)->GetAllRegions(allregions);


}

void KbGroup::Load(soci::session& sociDB) {
	sociDB << "SELECT group_name, group_desc FROM groups WHERE group_id=" <<DbID(), into(name), into(desc);
}



void KbGroup::AssociateGene(KbRegion* region) {
	if (isDiseaseDependent)
		region->InsertDD(groupType);
	else
		region->InsertDI(groupType);
	regions.insert(region);
}
uint KbGroup::AssociateGenes(soci::session& sociDB, map<uint, KbRegion*>& genes, SnpManager& snps, map<uint, string>& aliases, uint popID) {
	if (genesAssociated)
		return 0;

	cerr<<"We shouldn't ever associate genes through this function\n";
	assert(0);

	uint geneCount = 0;
	//We have to try to do this for children as well
	set<KbGroup*>::iterator child = children.begin();
	set<KbGroup*>::iterator childrenEnd = children.end();
	while (child != childrenEnd)
		geneCount += ((KbGroup*)*child++)->AssociateGenes(sociDB, genes, snps, aliases);


	map<uint, string>::iterator notFound = aliases.end();

	//Eventually, this will be updated to reflect the user's population
	//rowset<row> rs = (sociDB.prepare << "SELECT regions.gene_id, ensembl_id, chrom, start, end, description FROM (SELECT * FROM region_bounds WHERE population_id=:popID) NATURAL JOIN regions INNER JOIN group_associations ON (regions.gene_id=group_associations.gene_id) WHERE group_id=:id", use(popID), use(DbID()));
	rowset<row> rs = (sociDB.prepare << "SELECT gene_id FROM  group_associations  WHERE group_id=:id", use(popID), use(DbID()));
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint gene_id = row.get<int>(0);
		geneCount++;
		KbRegion *region = NULL;
		if (genes.find(gene_id) == genes.end()) {
			assert(0);
			uint start = row.get<int>(3);
			uint stop = row.get<int>(4);
			string chrom = row.get<string>(2);
			string ensembl = row.get<string>(1);
			string desc = "";		//row.get<string>(5);
			region = new KbRegion(gene_id, start, stop, chrom.c_str(), ensembl.c_str(),desc.c_str(), &snps);
			if (aliases.find(gene_id) != notFound)
				region->SetAlias(aliases[gene_id].c_str());
			region->AssociateSNPs();
			genes[gene_id] = (KbRegion*)region;
		}
		if (region == NULL)
			region = genes[gene_id];
		AssociateGene(region);
	}

	genesAssociated = true;

	return geneCount;
}

uint KbGroup::AssociateSNPs(SnpManager& snps) {
	set<KbGroup*>::iterator itr = children.begin();
	set<KbGroup*>::iterator end = children.end();

	uint snpCount = 0;
	while (itr != end) {
		snpCount += (*itr)->AssociateSNPs(snps);
		itr++;
	}
	return snpCount;
}


void KbGroup::MarkProcessed(bool isProcessed) {
	processed = isProcessed;
		
	set<KbGroup*>::iterator itr = children.begin();
	set<KbGroup*>::iterator end = children.end();
	while (itr != end)
		(*itr++)->MarkAsProcessed(isProcessed);

	set<KbRegion*>::iterator rItr = regions.begin();
	set<KbRegion*>::iterator rEnd = regions.end();

	while (rItr != rEnd)
		(*rItr++)->MarkAsProcessed(isProcessed);
}


void KbGroup::LoadAssociations(soci::session& sociDB, map<uint, KbGroup*>& groups, set<uint> parents) {
	if (processed)
		return;
	rowset<row> rs = (sociDB.prepare << "SELECT child_id FROM group_relationships WHERE parent_id=:id", use(DbID()));

	//We want to update local parents with those inside this calling function
	//AddParents(parents);

	//We want to add ourself to children parent lists
	//parents.insert(id);
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint groupID = row.get<int>(0);
		KbGroup *newGroup = NULL;

		if (groupID != DbID()) {
			if (groups.find(groupID) == groups.end()) {
				newGroup = new KbGroup(groupID, groupType, isDiseaseDependent, "", "");
				((KbGroup*)newGroup)->Load(sociDB);
				groups[groupID] = newGroup;
			}
			if (newGroup == NULL)
				newGroup = groups[groupID];
			children.insert(newGroup);
			newGroup->LoadAssociations(sociDB, groups, parents);
		}
	}
	MarkAsProcessed(true);
}

uint KbGroup::GraphAssociations(ostream& os, uint maxGeneCount) {
	static std::map<uint,uint> colors;
	stringstream ss;

	uint snpCount = 0;

	if (CollapseAssociationReport) {
		set<KbRegion*> allRegions;
		GetAllRegions(allRegions);
		if (allRegions.size() <= maxGeneCount) {
			//ss<<CommonName()<<" +("<<allGroups.size()<<", "<<allRegions.size()<<")"<<"\n";
			set<KbRegion*>::iterator rItr = allRegions.begin();
			set<KbRegion*>::iterator rEnd = allRegions.end();
			uint prevID = 0;//, firstID=0;
			while (rItr!=rEnd) {
				KbRegion* region = *rItr++;
				if (region->DdGroups().size() > 0) 
					snpCount++;
				ss<<"\t"<<region->DbID()<<" [label=\""<<region->CommonName()<<"\",shape=ellipse,color=pink];\n"
						<<"\t"<<DbID()<<"->"<<region->DbID()<<"\n";
				prevID = region->DbID();
			}
			if (snpCount > 0) {
				os<<"\t"<<DbID()<<" [label=\""<<CommonName()<<"\",shape=square,color=green];\n"<<ss.str();
			}
			return snpCount;
		}

	}

	set<KbGroup*>::iterator itr = children.begin();
	set<KbGroup*>::iterator end = children.end();
	stringstream childStream;
	while (itr != end) {
		uint snps = (*itr)->GraphAssociations(childStream, maxGeneCount);
		if (snps > 0)
			ss<<"\t"<<DbID()<<"->"<<(*itr)->DbID()<<";\n";
		snpCount+=snps;
		itr++;
	}

	uint prevID = 0;//, firstID=0;
	if (!CollapseAssociationReport || (CollapseAssociationReport && regions.size() < maxGeneCount)) {

		set<KbRegion*>::iterator rItr = regions.begin();
		set<KbRegion*>::iterator rEnd = regions.end();

		uint unIDdGenes = 0;
		//stringstream ss;
		while (rItr!=rEnd) {
			KbRegion* region = *rItr++;
			bool isDD=region->DdGroups().size() > 0;
			if (isDD)
				snpCount++;
			else
				unIDdGenes++;

			uint col = 1;
			if (isDD) {
				if (colors.find(region->DbID()) == colors.end())
					colors[region->DbID()] = (colors.size()+2) % 12;
				col = colors[region->DbID()];
			}
			if (isDD || unIDdGenes < 5)
				ss<<"\t"<<DbID()<<region->DbID()<<" [label=\""<<region->CommonName()<<"\",shape=ellipse,fillcolor="<<col<<"];\n"
				<<"\t"<<DbID()<<"->"<<DbID()<<region->DbID()<<";\n";
			prevID = region->DbID();
		}
		if (unIDdGenes > 5)
			ss<<"\t"<<"zzz"<<DbID()<<" [label=\"... ("<<unIDdGenes-5<<")\",shape=ellipse,fillcolor=1];\n"
				<<"\t"<<DbID()<<"->zzz"<<DbID()<<";\n";
	}

	if (snpCount > 0) {
		//cout<<"--> "<<DbID()<<" label: "<<CommonName()<<"\t"<<ss.str()<<"\n";
		os<<"\t"<<DbID()<<" [label=\""<<CommonName()<<"\",shape=box,fillcolor=2];\n"<<ss.str();
	}

	return snpCount;

}
	uint KbGroup::ListAssociations(int tabCount, ostream& os, uint maxGeneCount) {
		stringstream ss;

		uint snpCount = 0;
		PrintTabs(tabCount, ss);

		if (CollapseAssociationReport) {
			set<KbRegion*> allRegions;
			GetAllRegions(allRegions);
			set<KbGroup*>  allGroups;
			GetAllGroups(allGroups);
			if (allRegions.size() <= maxGeneCount) {
				ss<<CommonName()<<" +("<<allGroups.size()<<", "<<allRegions.size()<<")"<<"\n";
				set<KbGroup*> allGroups;
				set<KbRegion*>::iterator rItr = allRegions.begin();
				set<KbRegion*>::iterator rEnd = allRegions.end();

				while (rItr!=rEnd)
					snpCount += (*rItr++)->ListAssociations(tabCount+1, ss);


				if (snpCount > 0)
					os<<ss.str();
				return snpCount;
			}

		}

		set<KbGroup*>::iterator itr = children.begin();
		set<KbGroup*>::iterator end = children.end();
		
		stringstream childStream;
		while (itr != end) {
			snpCount+=(*itr)->ListAssociations(tabCount+1, childStream, maxGeneCount);
			itr++;
		}

		if (!CollapseAssociationReport || (CollapseAssociationReport && regions.size() < maxGeneCount)) {

			set<KbRegion*>::iterator rItr = regions.begin();
			set<KbRegion*>::iterator rEnd = regions.end();

			while (rItr!=rEnd)
				snpCount += (*rItr++)->ListAssociations(tabCount+1, childStream);
		}

		if (snpCount > 0)
			os<<ss.str()<<CommonName()<<" -("<<children.size()<<", "<<snpCount<<")"<<"\n"<<childStream.str();
		
		return snpCount;
	
	}
}
}
