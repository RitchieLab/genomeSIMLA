/* 
 * File:   pathway.cpp
 * Author: torstees
 * 
 * Created on January 7, 2010, 3:22 PM
 */

#include "pathway.h"
#include <vector>
#include <iomanip>
#include "chromosome.h"
#include "utility/stats.h"
namespace Paris {

using namespace std;

bool Pathway::countFeatureRepeats = true;

Pathway::Pathway() : groupID(0), name(""), desc("") { }
Pathway::Pathway(uint groupID, const char *name, const char *desc) : groupID(groupID), name(name), desc(desc) { }
Pathway::Pathway(const Pathway& other) : groupID(other.groupID), name(other.name), desc(other.desc) { }
Pathway::~Pathway() { }

set<Gene*> Pathway::GetGenes() {
	return genes;
}


uint Pathway::GroupID() {
	return groupID;
}


std::string Pathway::Name() {
	return name;
}
std::string Pathway::Description() {
	return desc;
}

void Pathway::AddGene(Gene* gene) {
	genes.insert(gene);
}

void Pathway::GenerateSnpReport(std::ostream& os, std::map<std::string, Chromosome*>& chromosomes, bool writeHeader) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	if (writeHeader)
		os<<"RS ID,Chromosome,Position (bp),P-Value,Gene,Feature Start,Feature End,Pathway Name\n";

	while (itr != end) {
		Gene *gene = *itr++;

		std::set<Feature*> features = gene->GetFeatures();

		std::set<Feature*>::iterator fItr = features.begin();
		std::set<Feature*>::iterator fEnd = features.end();
		Chromosome *c					= chromosomes[gene->_chromosome];

		while (fItr != fEnd) {
			Feature *f					= *fItr++;

			std::map<uint, float> pvalues = f->GetPValues();
			std::map<uint, float>::iterator pvItr = pvalues.begin();
			std::map<uint, float>::iterator pvEnd = pvalues.end();

			while (pvItr != pvEnd) {
				uint pos						= pvItr->first;
				float pval					= pvItr->second;
				pvItr++;
				uint rsid					= c->GetSNP(pos);
				os<<"rs"<<rsid<<","
					<<gene->_chromosome<<","
					<<pos<<","
					<<pval<<","
					<<gene->Alias()<<","
					<<f->_begin<<","
					<<f->_end<<","
					<<Name()<<"\n";
			}
		}
	}
}

void Pathway::ReportVitals(std::ostream& os) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	Accum<float> geneSizeCounter(0);
	Accum<int> geneContentsSizeCounter(0);
	Accum<float> avgLdBlockSize(0);


	std::set<uint> snps;
	std::set<Feature*> features;
	
	while (itr != end) {
		snps.empty();
		(*itr)->CollectSNPs(snps);
		features = (*itr)->GetFeatures();
		
		std::set<Feature*>::iterator fItr = features.begin();
		std::set<Feature*>::iterator fEnd = features.end();

		uint countedFeatures = 0;
		while (fItr != fEnd) {
			countedFeatures+=1;
			if ((*fItr)->FeatureSize() > 0)
				avgLdBlockSize.AddValue((*fItr)->FeatureSize());
			fItr++;
		}

		if (countedFeatures > 0) {
			geneContentsSizeCounter.AddValue((int)snps.size());
			geneSizeCounter.AddValue((*itr)->_end - (*itr)->_begin);
		}
		itr++;
	}

	os<<name<<"\t"<<groupID<<"\t"<<geneSizeCounter.GetMean()<<"\t"<<geneSizeCounter.GetStdDev()
				<<"\t"<<geneContentsSizeCounter.GetMean()<<"\t"<<geneContentsSizeCounter.GetStdDev()
				<<"\t"<<avgLdBlockSize.GetMean()<<"\t"<<avgLdBlockSize.GetStdDev()<<"\n";
}

void Pathway::ReportGenesAndFeatures(std::ostream& os) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	os<<"Gene List: "<<groupID<<"\t"<<name<<"\n";
	while (itr != end) {
		Gene* gene = *itr++;
		cerr<<"\t"<<gene->_chromosome<<"\t"<<gene->EnsemblID()<<"\t"<<gene->CountComplexFeatures()
			<<"\t"<<gene->CountSimpleFeatures()<<"\t"<<gene->CountEmptyFeatures()<<"\n";

	}
}

uint Pathway::TotalSNPs() {
	std::set<uint> snps;

	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	while (itr != end) {
		Gene *gene = *itr++;
		gene->CollectSNPs(snps);
	}
	return snps.size();
}
void Pathway::DetailedReport(std::map<std::string, Chromosome*>& chroms, const char *prefix, std::ostream& os) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();
	std::stringstream ss;
	ss<<prefix<<","<<name;
	while (itr != end) {
		Gene *gene = *itr++;
		std::map<uint, uint> snps = chroms[gene->_chromosome]->GetSNPs();
		gene->DetailedReport(snps, ss.str().c_str(), os);
	}
}

Analyzer::Result Pathway::RunAnalysis(std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose) {
	//ReportGenesAndFeatures(cerr);
	if (verbose) {
		cerr<<"Pathway: "<<name; cerr.flush();
	}

	Analyzer::Result result = analyzer.Analyze(groupID, genes, bins, permutationCount, verbose);

	return result;
}

Analyzer::Result Pathway::AnalyzeNegativeControl(std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose) {
	//ReportGenesAndFeatures(cerr);
	if (verbose) {
		cerr<<"Pathway: "<<name; cerr.flush();
	}

	Analyzer::Result result = analyzer.AnalyzeForNegativeControls(groupID, genes, bins, permutationCount, verbose);

	return result;
}



}

