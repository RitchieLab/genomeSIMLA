//
// C++ Interface: omnibusdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_REPORTINGOMNIBUSDISTRIBUTION_H
#define GENETICS_REPORTINGOMNIBUSDISTRIBUTION_H
#include "permutationtest.h"
namespace Genetics {

namespace Reporting {

/**
 * @brief Builds a test distribution based on the omnibus p-test
 * The omnibus test keeps only the best performing statistic for a given set of randomized statuses

 * @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class OmnibusDistribution : public PermutationTestDist {
public:
	~OmnibusDistribution()  {}
    OmnibusDistribution(uint seed, uint testCount, BasicLog *log = NULL) : PermutationTestDist(log), isSorted(false), seed(seed) {
		for (uint i=0; i<testCount; i++) 
			scores.push_back(TestEvaluation(i));
	}

	/**
 	 * @brief Provide a snp for consideration
	 */
	void Append(uint test, SnpAligned *snp);

	/**
	 * @brief Appends a new score to the distribution
	 */
	void Append(float score, uint test, SnpAligned *snp);

	/**
	 * @brief Return the pvalue associated with the current score in relationship to all other scores in the distribution
	 */
	float GetPValue(uint modelsize, float score);
	void ReportConfig(ostream *os);				///<Report on the configuration
	void Sort() ;								///<Perform sort

	bool isSorted;								///<Test to see if the results have been sorted

	void Report();								///<Generate a simple report

	FloatArray scores;							///<The array of scores 
	uint seed;									///<The seed associated with the set of runs
	void Append(float score, uint modelSize, uint test, const char *label);
};

inline
void OmnibusDistribution::Report() {
	ostream *os = NULL;
	if (log && log->GetStream())
		os = log->GetStream();

	if (os) {
		uint arrayCount=scores.size();
	
		for (uint h=0; h<arrayCount; h++) {
			if (scores[h].snp) {
				*os<<h<<" ("<<scores[h].snp->GetLabel()<<") "<<scores[h].score<<"\n";
			}
		}
	}
}

inline
void OmnibusDistribution::ReportConfig(ostream *os) {
	if (os == NULL)
		os = &cout;
	*os<<"Omnibus Permutation Test Distribution\n";
	*os<<"Random Seed: "<<seed<<endl;
	*os<<"Test Count: "<<scores.size()<<endl;
}


inline
void OmnibusDistribution::Append(float score, uint modelSize, uint test, const char *label) {
	assert(test < scores.size()); 
	isSorted=false;

	if (scores[test].score<score) 
		scores[test].Evaluate(label, score);

	isSorted=false;
}

inline
void OmnibusDistribution::Append(float score, uint test, SnpAligned *snp) {
	scores[test].Evaluate(snp, score);
	isSorted=false;
}

inline
void OmnibusDistribution::Append(uint test, SnpAligned *snp) {
	scores[test].Evaluate(snp);
	isSorted=false;		
}

inline
void OmnibusDistribution::Sort() {
	sort(scores.begin(), scores.end());

	isSorted=true;
}



}

}

#endif
