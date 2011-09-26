/* 
 * File:   analyzer.cpp
 * Author: torstees
 * 
 * Created on February 8, 2010, 10:30 AM
 */

#include "analyzer.h"
#include "feature.h"
#include <iostream>
#include <iomanip>
#include "assert.h"

using namespace std;

namespace Paris {

Utility::Random &Analyzer::rnd = Utility::Random::globalGenerator;
std::pair<float, float> Analyzer::refinementThreshold = std::pair<float, float>(0.04, 0.06);
Analyzer::Analyzer() {
}

Analyzer::Analyzer(const Analyzer& orig) {
}

Analyzer::~Analyzer() {
}

Analyzer::Result Analyzer::Analyze(uint reportID, std::set<Gene*>& genes, std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose) {
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	std::multiset<uint> binIDs;						///< This will tell us which bins to draw from
	std::set<uint> featureIDs;							///< Prevent from resampling one from the original source

	uint sigCount = 0;

	set<Feature*> usedFeatures;

	Analyzer::Result result(reportID);

	//Set up a pvalue of 1.0 to short circuit pathways with 0 significant features
	result.totPerms = permutationCount;
	result.sigPerms = permutationCount;
	result.geneCount = genes.size();

	while (itr != end) {
		Gene *gene = *itr;
		//binIDs will hold N copies of each of the bin IDs present, allowing us to
		//query how many of each bin we need to test
		gene->GetFeatureMap(binIDs, featureIDs);

		//Set up the number of features that have 1 or more signficant pvalues in them
		uint sigMembers = gene->CountSignificantMembers(usedFeatures, result.simpleFeatures, result.complexFeatures, result.sigSimple, result.sigComplex);
		if (sigMembers > 0)
			result.sigGenes++;

		sigCount += sigMembers;
		itr++;
	}

	if (verbose)
		cerr<<"\tGenes: "<<setw(5)<<genes.size()
			 <<"\tFeatures: "<<setw(5)<<binIDs.size()
			<<"\tSigificant: "<<setw(5)<<sigCount<<"..."; cerr.flush();

	//Short circuit to avoid testing situations that can never be less than 1.0
	if (sigCount > 0) {
		uint significantPTests=0;
		for (uint i=0; i<permutationCount; i++) {
			uint permutation = 0;
			//Draw without replacement-we don't want to allow these to be used in
			//permutations
			std::set<uint> untouchables = featureIDs;

			std::multiset<uint>::iterator binID = binIDs.begin();
			std::multiset<uint>::iterator binEnd = binIDs.end();

			while (binID != binEnd) {
				vector<Feature*>& features = bins[*binID];
				Feature *f = GetNextBin(untouchables, features);
				if (f->CountSignificantMembers() > 0)
					permutation++;
				untouchables.insert(f->id);
				binID++;
			}


			

			//Permutation is only significant if it has more sig. features than the original
			if (permutation > sigCount)
				significantPTests++;
		}
		result.sigPerms = significantPTests;
		if (verbose)
			cerr<<" Sig PTests: "<<setw(8)<<significantPTests<<" -> ";

	}
	else
		if (verbose)
		  cerr<<" Sig PTests: "<<setw(8)<<0.0<<" -> ";
	if (verbose)
		cerr<<result.PValue()<<"\n";
	return result;
}



Analyzer::Result Analyzer::AnalyzeForNegativeControls(uint reportID, std::set<Gene*>& genes, std::map<uint, std::vector<Feature*> >& bins, uint permutationCount, bool verbose) {
	verbose = false;
	std::set<Gene*>::iterator itr = genes.begin();
	std::set<Gene*>::iterator end = genes.end();

	std::multiset<uint> binIDs;						///< This will tell us which bins to draw from
	std::set<uint> featureIDs;							///< Prevent from resampling one from the original source

	uint sigCount = 0;

	set<Feature*> usedFeatures;

	Analyzer::Result result(reportID);

	//Set up a pvalue of 1.0 to short circuit pathways with 0 significant features
	result.totPerms = permutationCount;
	result.sigPerms = permutationCount;
	result.geneCount = genes.size();
	
	while (itr != end) {
		Gene *gene = *itr;
		//binIDs will hold N copies of each of the bin IDs present, allowing us to
		//query how many of each bin we need to test
		gene->GetFeatureMap(binIDs, featureIDs);

		//Set up the number of features that have 1 or more signficant pvalues in them
		uint sigMembers = gene->CountSignificantMembers(usedFeatures, result.simpleFeatures, result.complexFeatures, result.sigSimple, result.sigComplex);
		if (sigMembers > 0)
			result.sigGenes++;

		sigCount += sigMembers;
		itr++;
	}

	result.sigSimple = 0;
	result.sigComplex = 0;


	sigCount = 0;
	std::set<uint> visited = featureIDs;



	//Short circuit to avoid testing situations that can never be less than 1.0
	std::multiset<uint>::iterator binID = binIDs.begin();
	std::multiset<uint>::iterator binEnd = binIDs.end();

	std::set<uint> localFeatures;
	//Build up the local random pathway
	while (binID != binEnd) {
		vector<Feature*>& features = bins[*binID];
		Feature *f = GetNextBin(visited, features);
		if (f->CountSignificantMembers() > 0) {
			sigCount++;
			if (*binID == 0)
				result.sigSimple++;
			else
				result.sigComplex++;
		}
		localFeatures.insert(f->id);	//To avoid having any of the local features appear in permuted pathways
		visited.insert(f->id);			//To make sure we don't have any redundant features present
		binID++;
	}
	if (verbose)
		cerr<<"\tGenes: "<<setw(5)<<genes.size()<<"\tFeatures: "<<setw(5)<<binIDs.size()<<"\tSigificant: "<<setw(5)<<sigCount<<"..."; cerr.flush();

	uint significantPTests=0;
	for (uint p=0; p<permutationCount; p++) {
		uint permutation = 0;

		std::set<uint> untouchables = localFeatures;

		binID = binIDs.begin();
		binEnd = binIDs.end();
		while (binID != binEnd) {
			uint bIdx = *binID;
			vector<Feature*>& features = bins[bIdx];
			Feature *f = GetNextBin(untouchables, features);
			if (f->CountSignificantMembers() > 0)
				permutation++;
			untouchables.insert(f->id);
			binID++;
		}
		//Permutation is only significant if it has more sig. features than the original
		if (permutation > sigCount)
			significantPTests++;
	}
	cout<<"\t"<<significantPTests<<"\n";
	result.sigPerms = significantPTests;

	if (verbose)
		cerr<<result.PValue()<<"\n";
	return result;
}



Feature *Analyzer::GetNextBin(std::set<uint>& visited, std::vector<Feature*>& features) {
	Feature *f = 0;
	//assert(visited.size() < features.size());
	uint attempts = 1000000;
	do {
		uint fIdx = this->rnd((int)features.size());
		f = features[fIdx];
		attempts--;
		//If this fails, then we have run out of unique features
		assert(attempts > 0);
	} while (visited.find(f->id) != visited.end());
	return f;
}

}

