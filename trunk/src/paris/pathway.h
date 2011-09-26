/* 
 * File:   pathway.h
 * Author: torstees
 *
 * Created on January 7, 2010, 3:22 PM
 */

#ifndef _PATHWAY_H
#define	_PATHWAY_H
#include <iostream>
#include <set>
#include <string>
#include <assert.h>
#include <math.h>
#include "utility/random.h"
#include "utility/strings.h"
#include "analyzer.h"


namespace Paris {

class Chromosome;


class Pathway  {
public:
		Pathway();
		Pathway(uint groupID, const char *name, const char *desc);
		Pathway(const Pathway& other);
		~Pathway();

		uint GroupID();											///< Pathway ID (not meaningful out of the bio-settings database)
		std::string Name();										///< Name of the pathway
		std::string Description();								///< Description of the pathway

		uint TotalSNPs();
		void AddGene(Gene* gene);

		void ReportGenesAndFeatures(std::ostream& os);
		Analyzer::Result AnalyzeNegativeControl(std::map<uint, std::vector<Feature*> >& bins, uint permuationCount, bool verbose=true);
		Analyzer::Result RunAnalysis(std::map<uint, std::vector<Feature*> >& bins, uint permuationCount, bool verbose=true);
		void DetailedReport(std::map<std::string, Chromosome*>& chroms, const char *prefix, std::ostream& os);
		static bool countFeatureRepeats;						///< Used to turn on/off how we treat Features that appear more than once in a pathway

		std::set<Gene*> GetGenes();
		void GenerateSnpReport(std::ostream& os, std::map<std::string, Chromosome*>& chromosomes, bool writeHeader = false);
		void ReportVitals(std::ostream& os);
protected:
		uint groupID;												///< Pathway ID (not meaningful out of the bio-settings database)
		std::string name;											///< Name of the pathway
		std::string desc;											///< Description of the pathway
		std::set<Gene*> genes;									///< Genes associated with the pathway
		Analyzer analyzer;										///< This is responsible for performing the analyzis
};




}
#endif	/* _PATHWAY_H */

