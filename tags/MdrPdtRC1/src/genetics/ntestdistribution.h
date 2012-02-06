//
// C++ Interface: ntestdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONNTESTDISTRIBUTION_H
#define GENETICS_EVALUATIONNTESTDISTRIBUTION_H
#include "permutationtest.h"
namespace Genetics {

namespace Reporting {

/**
 * @brief The original approach for p-tests in which each loci gets it's own distribution 
   @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>

 */
class NTestDistribution : public PermutationTestDist  {
public:
	~NTestDistribution() {}
    NTestDistribution(uint seed, uint snpCount, uint testCount, BasicLog *log) : PermutationTestDist(log), isSorted(false), scores(snpCount), seed(seed) {
		srand(seed);
		for (uint h=0; h<snpCount; h++) 
			scores[h].insert(scores[h].begin(), testCount, TestEvaluation(0.0));
	}
	void ReportConfig(ostream *os);
	/**
 	 * Provide a snp for consideration
	 */
	void Append(uint test, SnpAligned *snp);
	void Append(float score, uint test, SnpAligned *snp);

	//Assumes that the vector has been sorted
	float GetPValue(uint modelSize, float score);

	void Sort() ;

	bool isSorted;
	vector <FloatArray> scores;	

	void Report();
	uint seed;
};



}

}

#endif
