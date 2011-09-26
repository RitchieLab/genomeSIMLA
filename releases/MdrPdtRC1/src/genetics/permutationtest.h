//
// C++ Interface: permutationtest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_REPORTINGPERMUTATIONTEST_H
#define GENETICS_REPORTINGPERMUTATIONTEST_H

#include "snprecipient.h"
#include "snpcontainer.h"
#include <iomanip>

namespace Genetics {
namespace Reporting {

/**
@brief Responsible for the compiling the permutation test report

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class TestEvaluation {
public:
	float score;
	SnpAligned *snp;
	
	TestEvaluation(float score, SnpAligned *snp) : score(score), snp(snp) {
		snp->IncrementInstanceCount();
	}
	TestEvaluation() : score(0.0), snp(NULL) {}	
	TestEvaluation(float score) : score(score), snp(NULL) {}

	TestEvaluation(const TestEvaluation& other) {
		score=other.score;
		snp=other.snp;
		if (snp)
			snp->IncrementInstanceCount();
	}
	
	TestEvaluation &operator=(const TestEvaluation &other) {
		score=other.score;
		snp=other.snp;
		if (snp)
			snp->IncrementInstanceCount();
		return *this;
	}		

	TestEvaluation &operator=(SnpAligned *snp) {
		this->snp=snp;
		snp->IncrementInstanceCount();
		return *this;
	}

	void Evaluate(SnpAligned *snp) {
		float newScore = snp->GetLastMdEval();
		if (newScore > score) {
			score = newScore;	
			if (this->snp)
				snp->ReduceInstanceCount();
			this->snp=snp;
			snp->IncrementInstanceCount();
			cout<<"Distribution receiving new snp: "<<snp->GetID()<<" "<<snp->GetLabel()<<" "<<score<<"\n";
		}
	}

	void Evaluate(SnpAligned *snp, float newScore) {
		if (newScore > score) {
			score = newScore;	
			if (this->snp)
				snp->ReduceInstanceCount();
			this->snp=snp;
			snp->IncrementInstanceCount();
		}
	}	

	bool operator<(const TestEvaluation& other) const {
		return score>other.score;
	}

	~TestEvaluation() {
		if (snp)
			snp->ReduceInstanceCount();
	}
};

typedef vector<TestEvaluation> FloatArray;
class PermutationTestDist {
public:
	PermutationTestDist(BasicLog *log = NULL) : curPos(0), log(log) {}
	virtual ~PermutationTestDist() {}
	virtual void ReportConfig(ostream *os) = 0;
	virtual void Report() = 0;
	/**
 	 * Provide a snp for consideration
	 */
	virtual void Append(uint test, SnpAligned *snp)=0;
	virtual void Append(float score, uint test, SnpAligned *snp)=0; 
	virtual void Append(float score, SnpAligned *snp) { Append(score, curPos++, snp); }

	//Assumes that the vector has been sorted
	virtual float GetPValue(uint modelSize, float score)=0;

	virtual void Sort()=0;
	uint curPos;				///<The number of appends thus far if each done in order
	BasicLog *log;				///<The log used to report the distrubtion details to

};


}

}

#endif
