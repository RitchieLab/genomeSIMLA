/* 
 * File:   knowledgetree.cpp
 * Author: torstees
 * 
 * Created on January 7, 2010, 9:49 AM
 */

#include "knowledgebase.h"
#include "pathway.h"
#include <string>
#include <iomanip>
#include <vector>
#include "magicdb.h"
#include "utility/stat.h"
using namespace soci;
using namespace std;

namespace Paris {

uint KnowledgeBase::AnalysisRepCount		= 1000;

KnowledgeBase::KnowledgeBase() : id(0) { }
KnowledgeBase::KnowledgeBase(uint id, const char *name) : id(id), name(name) { }
KnowledgeBase::KnowledgeBase(const KnowledgeBase& other) : id(other.id), name(other.name), pathways(other.pathways) { }
KnowledgeBase::~KnowledgeBase() { cerr<<"Goodbye knowledge!\n";}


void KnowledgeBase::Purge() {
	std::map<uint, Pathway*>::iterator itr = pathways.begin();
	std::map<uint, Pathway*>::iterator end = pathways.end();

	while (itr != end) 
		delete (itr++)->second;

}

Pathway *KnowledgeBase::GetPathway(uint id) {
	return pathways[id];
}
Pathway *KnowledgeBase::GetPathway(const char *pathway) {
	return pathwayLookup[pathway];
}

uint KnowledgeBase::GetID() {
	return id;
}

string KnowledgeBase::GetName() {
	return name;
}

uint KnowledgeBase::LoadKnowledge(soci::session& sociDB, std::map<std::string, Chromosome*>& chromosomes, ostream& os) {

	sociDB.prepare <<"SELECT group_type FROM group_type WHERE group_type_id=:id", use(id), into(name);
	rowset<row> rs = (sociDB.prepare << "SELECT group_id, group_name, group_desc, gene_id, chrom FROM groups NATURAL JOIN group_associations NATURAL JOIN regions WHERE group_type_id=:typeID ORDER BY chrom, group_id", use(id));
	uint count = 0;
	assert(chromosomes.size() > 0);
	Chromosome *chrom = chromosomes.begin()->second;

	Purge();

	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint groupID = row.get<int>(0);
		string name = row.get<string>(1);
		string desc = row.get<string>(2);
		uint geneID	= row.get<int>(3);
		string c = row.get<string>(4);
		if (chrom && chrom->ID() != c)
			chrom = NULL;
		if (chromosomes.find(c) != chromosomes.end())
			chrom = chromosomes[c];
		else {
			cerr<<"We are having trouble finding chromosome '"<<c<<"'\n";
			std::map<std::string, Chromosome*>::iterator cItr = chromosomes.begin();
			//while (cItr!=chromosomes.end())
			//	cerr<<"'"<<cItr++->first<<"' ";
			//cerr<<"\n";
		}
		if (chrom) {
			Gene *gene = chrom->GetGene(geneID);
			assert(gene);
			gene->AddGroup(id, groupID);

			if (pathways.find(groupID) == pathways.end()) {
				Pathway *pathway = new Pathway(groupID, name.c_str(), desc.c_str());
				pathways[groupID] = pathway;

				assert(pathwayLookup.find(name) == pathwayLookup.end());
				pathwayLookup[name] = pathway;
				if (ParisResults::resultsDB)
					ParisResults::db.sociDB<<"INSERT INTO pathways VALUES (:pathwayID, :kbID, :name, :desc)", use(groupID), use(id), use(name), use(desc);
			}
			pathways[groupID]->AddGene(gene);
			os<<id<<"\t"<<name<<"\t"<<desc<<"\t"<<gene->id<<"\t"<<gene->EnsemblID()<<"\t"<<gene->_chromosome<<"\t"<<gene->_begin<<"\t"<<gene->_end;
			os<<"\n";
			count++;
		}
	}
	return count;
}

void KnowledgeBase::DetailedReport(std::map<std::string, Chromosome*>& chroms, std::ostream& os) {
	std::map<uint, Pathway*>::iterator itr = pathways.begin();
	std::map<uint, Pathway*>::iterator end = pathways.end();

	while (itr != end) {
		itr++->second->DetailedReport(chroms, name.c_str(), os);
	}
}

void KnowledgeBase::RunPVRefinement(std::map<uint, vector<Feature*> >& bins, uint pCount, Pathway* pathway, Analyzer::Result& result) {
	Utility::Accum<int> acc(AnalysisRepCount);

	cerr<<"Running refinement for pathway: "<<pathway->Name()<<" pvalue="<<result.GetPValue(); cerr.flush();
	for (uint i=0; i<AnalysisRepCount; i++) {
		Analyzer::Result res = pathway->RunAnalysis(bins, pCount, false);

		acc.AddValue(res.sigPerms);
	}

	result.SetAsRefinement(acc.Mean(), acc.StdDev());
	cerr<<"->\t"<<result.GetPValue()<<"\n";

}


void KnowledgeBase::WritePathwayStatistics(ostream& os) {
	std::map<uint, Pathway*>::iterator itr = pathways.begin();
	std::map<uint, Pathway*>::iterator end = pathways.end();

	while (itr != end) {
		os<<name<<"\t";
		itr->second->ReportVitals(os);
		itr++;
	}
}

struct Refinements {
	Pathway* pathway;
	Analyzer::Result result;
	Refinements(Pathway* p, Analyzer::Result r) : pathway(p), result(r) { }
	bool operator<(const Refinements& other) const {
		return result<other.result;
	}
};

uint KnowledgeBase::AnalyzeNegativeControl(std::map<uint, vector<Feature*> >& bins, uint pCount, float significance, std::multiset<Analyzer::Result>& scores) {
	uint significantPathways = 0;
	std::map<uint, Pathway*>::iterator itr = pathways.begin();
	std::map<uint, Pathway*>::iterator end = pathways.end();


	std::vector<Refinements> toBeRefined;
	while (itr != end) {
		Pathway *pathway = (itr++)->second;
		Analyzer::Result result = pathway->AnalyzeNegativeControl(bins, pCount);
		scores.insert(result);
		//We want to note any pathways that meet or fall below our significance threshold
		if (result.PValue() <= significance)
			significantPathways++;
	}
	
	return significantPathways;
}

uint KnowledgeBase::RunAnalysis(std::map<uint, vector<Feature*> >& bins, uint pCount, float significance, std::multiset<Analyzer::Result>& scores) {
	uint significantPathways = 0;
	std::map<uint, Pathway*>::iterator itr = pathways.begin();
	std::map<uint, Pathway*>::iterator end = pathways.end();


	std::vector<Refinements> toBeRefined;
	while (itr != end) {
		Pathway *pathway = (itr++)->second;
		Analyzer::Result result = pathway->RunAnalysis(bins, pCount);

		if (result.IsBorderline())
			toBeRefined.push_back(Refinements(pathway, result));
		else
			scores.insert(result);

		//We want to note any pathways that meet or fall below our significance threshold
		if (result.PValue() <= significance)
			significantPathways++;
	}

	std::vector<Refinements>::iterator pItr = toBeRefined.begin();
	std::vector<Refinements>::iterator pEnd = toBeRefined.end();
	if (toBeRefined.size() > 0)
		cerr<<"Performing "<<toBeRefined.size()<<" refinements ("<<AnalysisRepCount<<" repititions each)\n";
	while (pItr != pEnd) {
		Refinements r = *pItr++;
		RunPVRefinement(bins, pCount, r.pathway, r.result);
		scores.insert(r.result);
	}

	return significantPathways;
}

void KnowledgeBase::ResultsHeader(std::ostream& os) {
	os<<"KB ID\tKB Name\tPathway ID\tPathway Name\tTotal SNP Count\tDescription\tP-Value\tGene Count\tGene Count*\tComplex Feature Count\tComplex Feature Count**\tSimple Feature Count\tSimple Feature Count**\n";
}

void KnowledgeBase::ReportResults(std::multiset<Analyzer::Result>& scores, ostream& os) {
	std::multiset<Analyzer::Result>::iterator itr = scores.begin();
	std::multiset<Analyzer::Result>::iterator end = scores.end();

	cerr<<setw(10)<<"KB ID"<<setw(8)<<"(Name)"<<setw(8)<<"Path. ID"<<setw(20)<<"(Name)"<<setw(8)<<"P-Value"<<setw(8)
			 <<"Gene"<<setw(8)<<"Gene*"<<setw(8)<<"#F > 1"<<setw(8)<<"#F > 1*"<<setw(8)<<"#F = 1"
			 <<setw(8)<<"#F = 1*"<<setw(10)<<"Std Dev.*"<<" Pathway Desc."<<"\n";
	while (itr != end) {
		Pathway *pathway = pathways[itr->groupID];
		os<<id<<"\t"<<name<<"\t"<<itr->groupID<<"\t"<<pathway->Name()<<"\t"<<pathway->TotalSNPs()<<"\t"<<pathway->Description()<<"\t"<<itr->GetPValue()<<"\t"
		 <<itr->geneCount<<"\t"<<itr->sigGenes<<"\t"<<itr->complexFeatures<<"\t"<<itr->sigComplex<<"\t"<<itr->simpleFeatures<<"\t"<<itr->sigSimple;

		if (itr->isRefinement) 
			os<<"\t"<<itr->stddev;
		else
			os<<"\t ";
		os<<"\n";

		if (itr->PValue() <= 0.05) {
			cerr<<setw(10)<<id<<setw(8)<<name<<setw(8)<<itr->groupID<<setw(20)<<pathway->Name()<<setw(8)<<itr->GetPValue()<<setw(8)
				<<itr->geneCount<<setw(8)<<itr->sigGenes<<setw(8)<<itr->complexFeatures<<setw(8)<<itr->sigComplex<<setw(8)<<itr->simpleFeatures
				<<setw(8)<<itr->sigSimple;
			if (itr->isRefinement) 
				cerr<<setw(10)<<itr->stddev;
			else
				cerr<<setw(10)<<" ";
			cerr<<" "<<pathway->Description()<<"\n";
		}
		itr++;
	}
}

void KnowledgeBase::ResultsFooter(std::ostream& os, float datasetSignificance, float pathwaySignificance) {
	os<<"* genes that contain one or more feature with a pvalue at or below "<<datasetSignificance<<".\n**pvalues of "<<pathwaySignificance<<" or less. \n";
}
	
}
