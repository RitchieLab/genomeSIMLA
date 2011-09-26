//
// C++ Implementation: modelreport
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "modelreport.h"

namespace Biofilter {

namespace Reporting {
/*
ModelReport::ModelReport(std::ostream& os) : os(os), db(NULL) { }

ModelReport::~ModelReport()  { }

void ModelReport::SetDB(soci::session& db) {
	this->db = &db;
}

void ModelReport::SnpLookup(multimap<uint, uint>& snpLookup) {
	this->snpLookup = snpLookup;
}

void ModelReport::AddModel(uint snp1, uint snp2) {
	Model model;
	model.AddLocus(snp1);
	model.AddLocus(snp2);
	models.insert(model);
	//For now, we'll assume that we only add models once...so, if a SNP shows up in the model 
	//list more than once, it will get overwritten in the owners structure
	owners[snp1] = ModelOwners();
	owners[snp2] = ModelOwners();
}

void ModelReport::GenerateReport() {
	std::set<Model>::iterator itr = models.begin();
	std::set<Model>::iterator end = models.end();
	os<<"SNPS\t\t\tGenes\tGroups\tGenes\tGroups\n";
	while (itr != end) {
		stringstream buffer;
		Model model = *itr++;
		std::vector<uint> loci = model.GetLoci();
		os<<"[ ";
		for (size_t i=0; i<loci.size(); i++) {
			uint locus = loci[i];
			os<<locus<<" ";
//			if (owners[locus].ReportMetaGroups(buffer))
//				buffer<<"\t";
			if (owners[locus].ReportRegions(buffer))
				buffer<<"\t";
			if (owners[locus].ReportGroups(buffer))
				buffer<<"\t";
		}
		os<<"]\t"<<buffer.str()<<"\n";
	}
}

bool ModelReport::Process(KbMetaGroup* meta) {
	KbMetaGroup::iterator itr = meta->begin();
	KbMetaGroup::iterator end = meta->end();

	bool isPresent = false;
	while (itr != end) {
		isPresent = Process(meta, ((Group*)itr->second)) || isPresent;	
		itr++;
	}
	
	return isPresent;
}

bool ModelReport::Process(KbMetaGroup *meta, KbGroup* group) {
	KbGroup::iterator itr =  group->begin();
	KbGroup::iterator end =  group->end();

	bool success = false;
	while (itr != end) {
		KbGroup* g = (KbGroup*)*itr;
		success = Process(meta, g) || success;
	}

	set<KbRegion*> regions = group->GetRegions();
	set<KbRegion*>::iterator rItr = regions.begin();
	set<KbRegion*>::iterator rEnd = regions.end();

	while (rItr != rEnd)
		success = Process(meta, group, *itr++) || success;
		
	return success;
}

bool ModelReport::Process(KbMetaGroup* meta, Group* group, Region* region) {
	std::map<uint, ModelOwners>::iterator itr = owners.begin();
	std::map<uint, ModelOwners>::iterator end = owners.end();
	multimap<uint, uint>::iterator lkupEnd = snpLookup.end();
	bool isPresent = false;
	while (itr != end) {
		uint rsID = itr->first;
		multimap<uint, uint>::iterator snpItr = snpLookup.lower_bound(rsID);
		multimap<uint, uint>::iterator snpEnd = snpLookup.upper_bound(rsID);
		while (snpItr != snpEnd) {
			uint pos = snpItr->second;
			if (region->IsPresent(pos)) {
				isPresent = true;
				itr->second.AddRegion(region);
				itr->second.AddGroup(group);
				itr->second.AddMetaGroup(meta);
			}
			snpItr++;
		}
		itr++;
	}
	return isPresent;
}


ModelReportHTML::ModelReportHTML(const char*filename) : ModelReport(cout) {
	file.open(filename);
}
ModelReportHTML::~ModelReportHTML() {
}


void ModelReportHTML::GenerateReport() {
	std::set<Model>::iterator itr = models.begin();
	std::set<Model>::iterator end = models.end();
	file<<"<HTML><HEAD><TITLE>Biofilter Model Report</TITLE></HEAD>\n";
	file<<"<BODY>\n<TABLE CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES>\n";
	file<<"<TR bgcolor='#F3EFE0'><TH>SNP 1</TH><TH>Gene</TH><TH>Groups</TH><TH>SNP 2</TH><TH>Gene</TH><TH>Groups</TH></TR>\n";
	while (itr != end) {
		stringstream buffer;
		Model model = *itr++;
		std::vector<uint> loci = model.GetLoci();
		file<<"<TR>";
		for (size_t i=0; i<loci.size(); i++) {
			uint locus = loci[i];
			file<<"<TD>"<<HTML::LinkSnpReference(locus)<<"</TD>";
			file<<"<TD>";
			owners[locus].ReportRegionsWithLink(file);
			file<<"</TD><TD>";
			owners[locus].ReportGroupsWithLink(file);
			file<<"</TD>";
		}
		file<<"</TR>\n";
	}
	file<<"</TABLE></HTML>\n";
}

*/
}

}
