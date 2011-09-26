//
// C++ Interface: bioapplication
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERBIOAPPLICATION_H
#define BIOFILTERBIOAPPLICATION_H

#include "kbmetagroup.h"
#include "kbregion.h"
#include <set>
#include <iomanip>
#include "modelreport.h"
#include "snptogenemanager.h"
namespace Biofilter {

namespace Reporting {

using namespace std;

struct RegionReport {
	std::ostream& os;
	RegionReport(std::ostream& os) : os(os) { }
	void PrintHeader() {
		os<<setw(30)<<"Region Name"<<setw(8)<<"Chrom."
			<<setw(12)<<"BP Start"<<setw(12)<<"BP Stop"<<"\n";
	}
	void operator()(Knowledge::KbRegion* object) {
		os<<setw(30)<<object->CommonName()
		<<setw(8)<<object->Chromosome()
		<<setw(12)<<object->Start()
		<<setw(12)<<object->End()<<"\n";
	}
};


struct SnpToGeneMapping {
	std::ostream& os;
	std::ostream& osFailed;
	bool htmlFormat;
	SnpToGeneMapping(std::ostream& os, std::ostream& osFailed, bool htmlFormat = false)
		: os(os), osFailed(osFailed), htmlFormat(htmlFormat) { }
	void PrintHeader() {
		if (htmlFormat) {
			os<<"<HTML><HEAD><TITLE>Biofilter SNP to Gene Report</TITLE></HEAD>\n";
			os<<"<BODY><TABLE CELLSPACING=1 CELLPADDING=3 border=1 RULES=ROWS FRAME=HSIDES>\n";	
			os<<"\t<TR bgcolor='#F3EFE0'><TH>RS ID</TH><TH>Chrom.</TH><TH>Region Name</TH><TH>Region Aliases</TH></TR>\n";
		}
		else {
			os<<setw(12)<<"RS ID"<<setw(8)<<"Chrom."<<setw(30)<<"Region Name"<<"\tRegion Aliases\n";
		}
	}
	void Close() {
		if (htmlFormat)
			os<<"</TABLE></BODY></HTML>\n";
	}
 	void operator()(map<uint, Knowledge::KbRegion *>& regions, vector<uint>& snps, SnpManager& snpMgr) {
		//Build the Snp->Region map
		multimap<uint, Knowledge::KbRegion *> snpToRegion;
		map<uint, Knowledge::KbRegion*>::iterator r = regions.begin();
		map<uint, Knowledge::KbRegion*>::iterator rEnd = regions.end();
		while (r != rEnd) {
			set<SNP_Details> localSNPs;
			r->second->CollectSnpDetails(localSNPs);
			set<SNP_Details>::iterator snp = localSNPs.begin();
			set<SNP_Details>::iterator snpEnd = localSNPs.end();
			while (snp != snpEnd) {
				snpToRegion.insert(pair<uint, Knowledge::KbRegion*>(snp->rsID, r->second));
				snp++;
			}
			r++;
		}
		PrintHeader();
		//Print out the matches
		vector<uint>::iterator snp = snps.begin();
		vector<uint>::iterator sEnd = snps.end();
		while (snp != sEnd) {
			uint rsNumber = *snp;
			multimap<uint, Knowledge::KbRegion*>::iterator region = snpToRegion.lower_bound(rsNumber);
			multimap<uint, Knowledge::KbRegion*>::iterator regEnd = snpToRegion.upper_bound(rsNumber);
			if (region != regEnd) {
				if (htmlFormat) 
					os<<"\t<TR bgcolor='#FFFFFF'><TD ROWSPAN="<<snpToRegion.count(rsNumber)<<"><A HREF='http://www.ensembl.org/Homo_sapiens/Variation/Summary?source=dbSNP;v=rs"<<rsNumber<<"'>"
						<<rsNumber<<"</A></TD>";
				else
					os<<"rs"<<rsNumber;
				int count = 0;
				while (region!=regEnd) {
					if (htmlFormat) {
						if (count++ > 0)
							os<<"</TR><TR>";
						os<<"<TD>"<<region->second->Chromosome()<<"</TD><TD>"
							<<"<A HREF='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g="<<region->second->Name()<<"'>"
							<<region->second->Name()<<"</A></TD>";
					}
					else
						os<<"\t"<<region->second->Chromosome()<<"\t"<<region->second->CommonName()<<"\t";
					set<string> aliases = region->second->Aliases();
					set<string>::iterator aItr = aliases.begin();
					set<string>::iterator aEnd = aliases.end();
					if (htmlFormat) 
						os<<"<TD>";
					int aliasCount = 0;
					while (aItr != aEnd) {
						if (aliasCount++ > 0)
							os<<",";
						if (htmlFormat)
							os<<*aItr++;
						else  
							os<<*aItr++;
					}
					if (htmlFormat)
						os<<"</TD>";
					else
						os<<"\thttp://www.ensembl.org/Homo_sapiens/Gene/Summary?g="<<region->second->CommonName();
					region++;
				}
				if (htmlFormat)
					os<<"</TR>\n";
				os<<"\n";
			}
			else 
				osFailed<<"rs"<<rsNumber<<"\n";
			snp++;
		}
		Close();
	}
};
}
/**
@brief Performs various actions associated with biofilter. This class will not be dependant on a specific configuration method and should be useful to expose necessary functionality to any application that might use biofilter functions

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class BioApplication : public SnpManager {
public:
    BioApplication(const char *prefix = "", bool htmlReports = false);
    ~BioApplication();

	/**
	 * @Brief Initializes biofilter groups
	 * SNPs should be loaded into memory by this point, so the genes will collect
	 * their SNPs immediately. Empty genes and SNPs which aren't associated with
	 * one or more active groups will be removed from memory.
	 * @param inclusions Vector of Group/Type IDs which will be included. If empty, all groups and types will be included
	 * @param maxSizeForActive This sets max group size for inclusion. Groups will only be included if their size is below this value
	 * @param dbFilename This is the source database which contains group and variation IDs. 
	 */
	void InitBiofilter(const char *dbFilename, bool reportVersion=true);

	/**
	 * @brief Loads user's disease specific data into the group/region datastructures.
	 * This can be called more than once. A metagroup ID will be generated by the application during loading
	 * @param filename The source of the group/region information
	 * @param maxSizeForActive Maximum size (genes) before we shut off the group for model generation
	 * @return metaGroupID being used. I have no idea why they would want this 
	 * 
	 */
	//uint LoadUserDefinedGroupData(const char *filename, int maxSizeForActive);

	void ListPresentAssociations(uint maxGeneCount);

	void GraphPresentAssociations(const char *filename, uint maxGeneCount);

	/**
	 * @Brief details the coverage associated with a given set of SNPs and a set of genes
	 * @param geneList vector of genes to be initialized 
	 * @param os Stream to be written to
	 * This function assumes that the SNPs have already been loaded
	 */
	void DetailCoverage(std::vector<std::string>& genelist, std::vector<std::string>& snpFiles, bool detailedCovered);


	/**
	 * @Brief Initializes indexes in the database to allow for faster read. 
	 * @note Users should turn "Optimization" off when performing any sort of insertion, like loading the database with new group/LD imformation
	 */
	void PerformOptimization();
	
	/**
 	 * @Brief Drops indexes.
	 * @note this allows for faster insertions/deletions
	 */
	void StripOptimization();

	/**
	 * @brief Adds user defined grouping datasource
	 * @note This must be done prior to loading group data, since the contents of the files will be loaded then
	 */
	void AddUserDefinedGroup(const char *filename);


	/**
	 * @brief Produce models based on local group/snp information and store them into repo
	 * @param repo The destination for the snp-snp models
	 * @param os Stream where the system will write information
	 * @param maxGeneCount This is the maximum number of genes that a group can have in order to still produce models. 
	 * @note If a group exceeds maxGeneCount, the group will ask it's children to produce models, but not do it directly
	 */
	void ProduceModels(GeneGeneModelArchive& repo, std::ostream& os, int maxGeneCount);

	/**
	 * @Brief Generates a model summary report
	 */
	void SummarizeModelCounts(int maxGeneCount);
	/** 
	 * @brief Selectively loads genes into memory (snps must have already been loaded for snp association to work here)
	 * This function does not perform group associations.
	 * @param geneList Comma separated list of gene aliases
	 * @param pop Population ID associated with LD based region expansion
	 * @return number of genes that were matched from the list
	 */
	uint LoadRegions(const std::string& geneList, const char *pop);

	/**
	 * @Brief Loads group information into memory
	 * @param maxSizeForActive 
	 * @param includedGroups list of groups to be used. If this is empty, then all groups are retained
	 * @param pop Population ID associated with LD based region expansion
	 * @param prefGeneNames Filename with a list of region aliases. Any regions will be reported with these aliases, if present in the list
	 */
	void LoadGroupData(int maxSizeForActive, std::vector<uint>& includedGroups, std::ostream& os, const char *pop, const char *prefGeneNames = NULL);

	/**
	 * @Brief Builds a lookup table for geneID -> gene alias based on the contents of filename
	 * @param filename This is just a file with preferred gene aliases on separate lines
	 */
	std::map<uint, std::string> LoadRegionAlias(const char *filename);

	/**
	 * @Brief Loads LD according to the contents of file, ldConfiguration
	 * @param ldConfiguration This file specifies which population we are describing, the parameters associated with LD-Spline and the sources for the LD (along with the chromosome each file is associated with)
	 * @param variationFilename This is just the binary SNP information
	 */
	void ImportLD(const char *ldConfiguration, const char *variationFilename);
	void ImportLdSplines(const char *ldConfiguration);
	/**
	 * @brief Load merged data into rsConversion map, which is will be considered when loading SNP data from the file
	 */
	void LoadConversion();
	/**
	 * @BRief Returns the population ID from the database for a given string
	 */
	int GetPopID(const char *pop);
	std::string GetPopulationDesc(const char *pop);
	/**
	 * @Brief Produces a list of all population IDs in the database. These can then be used to associate LD with model generation
	 */
	void ListPopulationIDs();

	/**
	 * @Brief Produces a list of group IDs. These IDs can be used to restrict the search space used for model generation
	 */
	void ListGroupIDs(std::vector<std::string>& searchCriteria);

	void SetVariationFilename(const char *filename);

	/**
	 * @Brief Opens a model file and produces a text report (RS numbers, and portion of the implication index)
	 */
	//uint ReportOnModels(const char *modelFilename);

	/**
	 * @Brief Produce a report of all regions and groups/metagroups which directly contain one or more SNPs in a model
	 */
	//void RunReport(Reporting::ModelReport* report);
	/**
	 * @Brief Returns a pointer to the region object at geneID
	 * @param geneID The actual key from the database associated with a given region
	 */
	Knowledge::KbRegion *GetRegion(uint geneID);

	/**
	 * @Brief This is probably dead code...I think we are abandoning the idea of cross-group modelss
	 */
	static bool CrossGroupModelGen;				///< This is used to toggle model generation

	/**
	 * @Brief Displays the contents of a repository to stream, os.
	 * @param repo This contains the models to be displayed
	 * @param os Where we are writing the report
	 */
	//void DisplayContents(GeneGeneModelArchive& repo, std::ostream& os);

	/**
	 * @Brief Produce an overview of the SNPs and which regions they were found in
	 * @note Eventually, we'll want to add group information too.
	 */
	void SnpReport(std::ostream& os, std::ostream& failedSnps, std::vector<uint>& snps, bool writeHTML);

	/**
	 * @Brief Loads all region aliases associated with EntrezGene into the regions
	 */
	void LoadRegionAliases(const std::string& geneList = "ALL");
	void GetAliasMap(std::map<std::string, uint>& idLookup);
	void SetReportPrefix(const char *prefix);
	void UseHtmlReports(bool doUse);
	void ListMetaGroups(std::ostream& os);
	std::string GetReportLog();
	void ListGenes(std::ostream& os);
	static uint geneExtension;							///< Constant extension to apply where no LD is is use
	void ReportSNPsInGenes(std::vector<std::string>& genes, std::ostream& os);
	void BuildSnpToGene(SnpToGeneManager& snpsToGenes, std::vector<std::string>& geneList);
	uint LoadVariations(std::set<uint>& snps, const char *mergedReport, bool doReport);

	static bool PurgeOldIDs;							///< When true, we "rename" snps (and delete them). Otherwise, we let them both co-exist
protected:
	void InitGeneLookup(std::map<uint, Knowledge::KbRegion*>& regions);
	/**
	 * @Brief Loads regions into memory, based on comma separated file of gene aliases
	 * @param geneList The list of aliases of interest
	 * @Param regions This is where we store the regions (indexed by geneID)
	 * @Param os The stream to write report information to
	 * @param pop The character identifier associated with a given population
	 */
	uint LoadRegions(const std::string& geneList, std::map<uint, Knowledge::KbRegion*>& regions, std::ostream& os, const char *pop);

	/**
	 * @Brief Loads region information based on the population described by the text id, pop
	 */
 	uint LoadRegions(const char *pop);

	/**
	 * @brief Override current region boundaries with those found at population (pop)
	 * This should be used only to reassign. The population should be used during initial loading
	 */
	void SetRegionBoundaries(const char *pop, std::map<uint, Knowledge::KbRegion*>& regions);


	/**
	 * @Brief Loads disease dependent information into groups
	 * @param groupID The group ID to be associated wtih the new disease dependent group
	 * @param groupTypeID Group types are used for model generation..
	 * @param os To be written to during the load. This is for general reporting
	 * @param filename The file containing the disease dependant group information
	 * @param pop The id associated with the population ID (used for LD adjustment to gene boundaries)
	 */
	void LoadDiseaseDependent(uint &groupID, uint &groupTypeID, std::ostream& os, const char *filename, const char *pop);

	soci::session sociDB;					///< The database 
	std::map<uint, Knowledge::KbRegion*> regions;///< The genes that are associated with groups in memory
	std::map<uint, Knowledge::KbGroup*> groups;	///< Single lookup for each group object
	std::map<uint, Knowledge::KbMetaGroup*> metagroups;		///< The groups (PFam and KEGG are metagroups)


	/** Files where disease dependent data can be loaded */
	std::vector<std::string> diseaseDependentFiles;
	
	//These are used to produce IDs on the fly for stuff that isn't in our database
	uint maxGroupID;										///< highest group ID seen so far
	uint maxGroupTypeID;									///< highest group type ID seen so far
	uint maxRegionID;										///< highest region ID we've seen so far
	std::string reportPrefix;							///< Report Prefix
	bool htmlReports;										///< Indicate true for a preference for HTML reports
	std::string populationDesc;						///< Description relating to population used to LoadRegion
	std::stringstream reportLog;						///< Captures report filenames for announcing what was written where
	std::string dbFilename;								///< The database we are using
};




}

#endif
