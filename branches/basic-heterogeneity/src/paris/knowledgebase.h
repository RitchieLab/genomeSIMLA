/* 
 * File:   knowledgetree.h
 * Author: torstees
 *
 * Created on January 7, 2010, 9:49 AM
 */

#ifndef PARIS_KNOWLEDGETREE_H
#define	PARIS_KNOWLEDGETREE_H

#include <map>
#include <string>
#include <iostream>
#include "chromosome.h"
#include "pathway.h"

namespace Paris {


class KnowledgeBase {
public:
	KnowledgeBase();
	KnowledgeBase(uint id, const char *name);
	KnowledgeBase(const KnowledgeBase& orig);
	virtual ~KnowledgeBase();

	/**
	 * @brief Load up the data from the database, and grab associated genes from the chromsomes
	 */
	uint LoadKnowledge(soci::session& sociDB, std::map<std::string, Chromosome*>& chromosomes, std::ostream& os);

	/**
	 * @brief Run analysis for entire project
	 * @return Number of significant pathways
	 */
	uint RunAnalysis(std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, float significance, std::multiset<Analyzer::Result>& scores);
	/**
	 * @brief Run analysis for entire project
	 * @return Number of significant pathways
	 */
	uint AnalyzeNegativeControl(std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, float significance, std::multiset<Analyzer::Result>& scores);

	void RunPVRefinement(std::map<uint, std::vector<Feature*> >& bins, uint pCount, Pathway* pathway, Analyzer::Result& result);
	void DetailedReport(std::map<std::string, Chromosome*>& chroms, std::ostream& os);
	/**
	 * @brief Run analysis for a single group
	 * @return Number of permutations that exceeded the local score
	 */
	//uint RunAnalysis(uint groupID, std::map<uint, Feature*>& bins, uint permuationCount);

	void ReportResults(std::multiset<Analyzer::Result>& scores, std::ostream& os);
	void ResultsHeader(std::ostream &os);
	void ResultsFooter(std::ostream &os, float pvSig, float pathSig);
	Pathway *GetPathway(const char *pathwayName);
	Pathway *GetPathway(uint id);
	void WritePathwayStatistics(std::ostream& os);
	void Purge();

	uint GetID();

	std::string GetName();

	struct Score {
		Score() : groupID(0), pvalue(0.0) { }
		Score(uint groupID, float pvalue) : groupID(groupID), pvalue(pvalue) { }
		~Score() { }

		bool operator<(const Score& other) {
			return pvalue < other.pvalue;
		}

		uint groupID;
		float pvalue;
	};


	static uint AnalysisRepCount;
private:
	uint id;													///< Unique ID
	std::string name;										///< Name associated with the knowledge base (i.e. KEGG)
	std::map<uint, Pathway *> pathways;				///< List of group IDs found in genes
	std::map<std::string, Pathway*> pathwayLookup;		///< Lookup the pathway based on the pathway name
};




}
#endif	/* _KNOWLEDGETREE_H */

