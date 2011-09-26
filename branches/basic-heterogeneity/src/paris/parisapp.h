/* 
 * File:   parisapp.h
 * Author: torstees
 *
 * Created on January 5, 2010, 1:29 PM
 */

#ifndef PARIS_PARISAPP_H
#define	PARIS_PARISAPP_H

#include "chromosome.h"
#include <string>
#include <map>
#include <vector>
#include "knowledgebase.h"
namespace Paris {

class ParisApp {
public:
	struct DataPoint {
		DataPoint() : rsID(0), chrom(""), pvalue(0.0) { }
		DataPoint(uint rsID, const char *chrom, float pvalue) : rsID(rsID), chrom(chrom), pvalue(pvalue){ }
		~DataPoint() {}

		uint rsID;
		std::string chrom;
		float pvalue;
	};

	ParisApp();

	ParisApp(const ParisApp& orig);

	virtual ~ParisApp();

	/**
	 * @Brief Initializes biofilter groups
	 * SNPs should be loaded into memory by this point, so the genes will collect
	 * their SNPs immediately. Empty genes and SNPs which aren't associated with
	 * one or more active groups will be removed from memory.
	 * @param inclusions Vector of Group/Type IDs which will be included. If empty, all groups and types will be included
	 * @param maxSizeForActive This sets max group size for inclusion. Groups will only be included if their size is below this value
	 * @param dbFilename This is the source database which contains group and variation IDs.
	 */
	void InitKnowledge(const char *dbFilename);
	
	void InitKB(const char* popID, uint geneExpansion, Utility::StringArray& groups);
	
	void InitData(std::vector<ParisApp::DataPoint>& data, uint binSize);
	
	void SetReportPrefix(const char *prefix);
	
	uint InitSNPs(std::multimap<std::string, uint>& allSNPs);
	
	//void CleanRSIDs(std::map<uint, std::string>& snpList, const char *rsCleanReportFilename);
	void InvestigatePathway(uint permutationCount, float dataSig, float pathSig, const char *pathwayName, bool showAllPathways);
	uint InitBins(uint binSize);
	std::set<Pathway*> GetPathways(std::multimap<uint, uint>& groups);

	void RunAnalysis(const char *datafile, uint binSize, uint permutationCount, float datasetSignificance, float pathSignificance);
	
	void RunAnalysis(uint permutationCount, float datasetSignificance, float pathSignificance);

	std::map<std::string, uint> AnalyzeNegativeControl(uint permutationCount, float datasetSignificance, float pathSignificance);

	void ListGroupIDs();

	void WriteBinReport(const char *name);

	void WriteKbVitalsReport(const char *vitalsReport);

	void ReportName(const char *name);
	std::string ReportName();
private:
	soci::session sociDB;									///< The database
	std::string reportPrefix;								///< Report Prefix
	std::string reportName;									///< Optional name for the report to be used inside certain reports
	std::map<std::string, Chromosome*> chromosomes;	///< The chromosomes will manage the features and genes
	std::map<uint, KnowledgeBase*> knowledge;			///< id -> KB*
	std::map<uint, std::vector<Feature*> > bins;		///< binID -> Feature*
	std::string dbFilename;									///< The database we are using

};
}
#endif	/* _PARISAPP_H */

