//
// C++ Interface: modelreport
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTER_REPORTINGMODELREPORT_H
#define BIOFILTER_REPORTINGMODELREPORT_H
#include <set>
#include <map>
#include <string>
#include "kbmetagroup.h"
#include "kbgroup.h"
#include "html-helpers.h"

namespace Biofilter {

namespace Reporting {



/**
 * @brief Store only the groups and regions that directly contain this snp
 */
struct ModelOwners {
	std::set<Knowledge::KbGroup*> groups;
	std::set<Knowledge::KbRegion*> regions;
	std::set<Knowledge::KbMetaGroup*> metaGroups;

	void AddRegion(Knowledge::KbRegion* region) {
		regions.insert(region);
	}

	void AddGroup(Knowledge::KbGroup* group) {
		groups.insert(group);
	}

	void AddMetaGroup(Knowledge::KbMetaGroup* meta) {
		metaGroups.insert(meta);
	}
	bool ReportRegionsWithLink(std::ostream& os, std::string sep=std::string("|")) {
		std::set<Knowledge::KbRegion*>::iterator itr = regions.begin();
		std::set<Knowledge::KbRegion*>::iterator end = regions.end();
		os<<"<TABLE>";
		int count = 0;
		while (itr != end ){
			os<<"<TR><TD>"<<HTML::LinkGeneReference((*itr)->Name().c_str(), (*itr)->CommonName().c_str())<<"</TD>";
			itr++;
		}	
		os<<"</TABLE>";
		return count > 0;
	}		
	bool ReportRegions(std::ostream& os, std::string sep=std::string("|")) {
		std::set<Knowledge::KbRegion*>::iterator itr = regions.begin();
		std::set<Knowledge::KbRegion*>::iterator end = regions.end();

		int count = 0;
		while (itr != end ){
			if (count++ > 0)
				os<<sep;
			os<<(*itr)->CommonName();
			itr++;
		}
		return count > 0;
	}
	bool ReportGroupsWithLink(std::ostream& os) {
		std::set<Knowledge::KbGroup*>::iterator itr = groups.begin();
		std::set<Knowledge::KbGroup*>::iterator end = groups.end();
		int count = 0;
		os<<"<TABLE>";
		while (itr != end ){
			os<<"<TR><TD>"<<(*itr)->CommonName()<<"</TD></TR>";
			itr++;
		}
		os<<"</TABLE>";
		return count>0;
	}
		
	bool ReportGroups(std::ostream& os, std::string sep=std::string("|")) {
		std::set<Knowledge::KbGroup*>::iterator itr = groups.begin();
		std::set<Knowledge::KbGroup*>::iterator end = groups.end();
		int count = 0;
		while (itr != end ){
			if (count++ > 0)
				os<<sep;
			os<<(*itr)->CommonName();
			itr++;
		}
		return count>0;
	}

	bool ReportMetaGroups(std::ostream& os) {
		std::set<Knowledge::KbMetaGroup*>::iterator itr = metaGroups.begin();
		std::set<Knowledge::KbMetaGroup*>::iterator end = metaGroups.end();

		int count = 0;
		while (itr != end ){
			if (count++ > 0)
				os<<"|";
			os<<(*itr)->Name();
			itr++;
		}
		return count>0;
	}
};


/**
	@Brief Generates report describing the genes each model can be associated with
	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
/*class ModelReport{
public:
    ModelReport(std::ostream& os);
  	virtual ~ModelReport();

	void SetDB(soci::session& db);
	void SnpLookup(std::multimap<uint, uint>& snpLookup);

	//Quick 2 SNP model insert
	void AddModel(uint snp1, uint snp2);
	virtual void GenerateReport();

	bool Process(Knowledge::KbMetaGroup* meta, Knowledge::KbGroup* group, Knowledge::KbRegion* region);
	bool Process(Knowledge::KbMetaGroup* meta, Knowledge::KbGroup* group);
	bool Process(Knowledge::KbMetaGroup* meta);

protected:
	std::ostream& os;								///< where we are writing to
	soci::session* db;								///< The Database, in case we want to do some additional queries
	std::map<uint, ModelOwners> owners;				///< The owners for each endpoint (SNP)
	std::set<Model> models;							///< The actual models we are looking for
	std::multimap<uint, uint> snpLookup;			///< This is required to map rs numbers to positions, which is how everything is stored
};

class ModelReportHTML : public ModelReport {
public:
	ModelReportHTML(const char*filename);
	virtual ~ModelReportHTML();

	void GenerateReport();

protected:
	ofstream file;
};



*/

}

}

#endif
