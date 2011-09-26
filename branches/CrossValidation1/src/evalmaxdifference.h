//
// C++ Interface: evalbalancedaccuracy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEEVALMAXDIFFERENCE_H
#define ESEEVALMAXDIFFERENCE_H

#include "genetics/snpevaluation.h"
#include "utility/utility.h"
#include "genetics/snprepostxtsorted.h"
namespace ESE {

using namespace Utility;
using namespace Genetics;
using namespace Genetics::Reporting;
/**
 * 	@brief This module is for evaluating a repository for maximum difference
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class EvalMaxDifference : public Genetics::Evaluation::SnpEvaluation {
public:
    class Results;
	
	/**
	 * @brief Construction
	 * @param eval Pointer to the actual evaluation method.
	 */
	EvalMaxDifference();
	EvalMaxDifference(uint maxDiff, uint loci) : lociCount(loci), threshold(maxDiff) {}

	~EvalMaxDifference();

	/**
	 * @brief Adds a new fold to the set. 
	 * Each fold represents a set of weighted statuses to be used in identifying the top model
	 */
	void AppendStatus(CaseControlStatus &fold);

	/**
	 * @brief Evaluate a model over the ptests and the status set. 
	 * This represents a single evaluation as opposed to a x-validation one	
	 */
	bool EvaluateModel(SnpAligned *snp, uint position);

	void ReportResults(ostream *os);

	/**
	 * This represents the leading score seen for this module
	 */
	float GetTopScore(SnpAligned *&snp);
protected:
	/**
	 * @brief Evaluate a given snp for a specified status 
	 */
	uint Evaluate( SnpAligned *snp, const CaseControlStatus &status);
	//uint Evaluate(SnpAligned *snp, const StatusContainer &status);

	vector <CaseControlStatus> testingSet;
	vector <Results *> results;
	uint lociCount;
	uint threshold;
	
	
public:
/**
 * @brief Results associated with a single BA run
 */
struct Results {
	Results(uint threshold) : report(0), topModel(NULL), bestScore(0), threshold(threshold) {}
	~Results() {
		topModel->ReduceInstanceCount();
	}

	/**
	 * @brief Determine if the model is better than any other seen so far
	 */
	bool Evaluate(SnpAligned *snp, uint score) {
		bool isOK = false;
		if (score >= threshold) {
			isOK=true;
			report.Append(snp);
			if (score > bestScore) {
				bestScore = score;
				if (topModel)
					topModel->ReduceInstanceCount();
				topModel = snp;
				topModel->IncrementInstanceCount();
			}
		}
		return isOK;
	}

	SnpReposTxtSorted report;	///<Used to produce a basic report of the the top N models found
			
	SnpAligned *topModel;		///<Used to record the topmost model
	uint bestScore;				///<The current high score
	uint threshold;				///<The threshold associated with what is considered interesting

};
};

// TODO I need to set this up to work with x-validation
inline
float EvalMaxDifference::GetTopScore(SnpAligned *&snp) {
	uint count=results.size();
	uint bestResult =0;
	for (uint i=0; i<count; i++) { 
		if (bestResult < results[i]->bestScore) {
			bestResult=results[i]->bestScore;
			snp=results[i]->topModel;
		}
	}
	return (float)bestResult;
}

inline
EvalMaxDifference::EvalMaxDifference() : lociCount(0), threshold(0) {
	
}

inline
EvalMaxDifference::~EvalMaxDifference() {
	uint count=results.size();
	for (uint i=0; i<count; i++) 
		delete results[i];
	results.clear();
}

inline
void EvalMaxDifference::AppendStatus( CaseControlStatus &fold) {
	testingSet.push_back(fold);
	Results *newResults = new Results(threshold);
	results.push_back(newResults);
}


inline
bool EvalMaxDifference::EvaluateModel(SnpAligned *snp, uint position) {
	uint count=results.size();
	uint score=0;
	bool exceedsThreshold = true;
	for (uint i=0; i<count; i++) {
		score=Evaluate(snp, testingSet[i]);
		exceedsThreshold = exceedsThreshold & results[i]->Evaluate(snp, score);
	}
	return exceedsThreshold;
}

inline
uint EvalMaxDifference::Evaluate( SnpAligned *snp, const CaseControlStatus &status) {
	int affCount = 0;				///<The number of bits masked to affected
	int unaffCount = 0;
	int indivCount = 0;				///<Used to store the total number of 1s within a given bitmask
	int difference = 0;				///<Used to hold the current difference value
	uint mdValue=0;

	BitSetType a;
	BitSetType u;
	BitSetType curGenotype;	
	uint gtCount = snp->CountGenotypes();

	for (uint i=1; i<gtCount; i++) {
		curGenotype=snp->GetGenotype(i);
		a = curGenotype & status.affected;
		u = curGenotype & status.unaffected;

		affCount = a.count();
		unaffCount = u.count();
		indivCount = curGenotype.count();
		
		//It might be worth checking to see if this is slower than inverting the mask and counting the
		//number of unaffected individuals
		difference = abs(affCount+affCount - indivCount);
		//difference = abs(affCount - unaffCount);
		//Below is the same as affected - unaffected
		mdValue += difference;				//Absolute value;
		//cout<<GetLabel()<<" : "<<idx++<<" - "<<difference<<" "<<t<<" - "<<(*i & ~mask)<<"\n";	

	}
	return mdValue;
}
}

#endif

