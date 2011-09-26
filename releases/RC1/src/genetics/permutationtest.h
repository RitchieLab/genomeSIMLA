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
	uint id;					///<The test ID
	float score;				///<The score of the test
	SnpAligned *snp;			///<The snp holding this score
	string label;				///<The ID associated with the snp
	
	TestEvaluation(uint id, float score, SnpAligned *snp) : id(id), score(score), snp(snp) {
		label = snp->GetLabel();
		snp->IncrementInstanceCount();
	}
	TestEvaluation() : id(0), score(0.0), snp(NULL), label("") {}	
	TestEvaluation(uint id) : id(id), score(0.0), snp(NULL), label("") {}
	TestEvaluation(uint testID, float score, const char *label) : id(testID), score(score), snp(NULL), label(label) {}

	TestEvaluation(const TestEvaluation& other) {
		score=other.score;
		snp=other.snp;
		label=other.label;
		id=other.id;
		if (snp)
			snp->IncrementInstanceCount();
	}
	
	TestEvaluation &operator=(const TestEvaluation &other) {
		score=other.score;
		snp=other.snp;
		label=other.label;
		id=other.id;
		if (snp)
			snp->IncrementInstanceCount();
		return *this;
	}		

	TestEvaluation &operator=(SnpAligned *snp) {
		score=snp->GetLastMdEval();
		label=snp->GetLabel();
		this->snp=snp;
		snp->IncrementInstanceCount();
		return *this;
	}

	void Evaluate(SnpAligned *snp) {
		float newScore = snp->GetLastMdEval();
		if (newScore > score) {
			score = newScore;	
			label=snp->GetLabel();
			if (this->snp)
				snp->ReduceInstanceCount();
			this->snp=snp;
			snp->IncrementInstanceCount();
			cout<<"Distribution receiving new snp: "<<snp->GetID()<<" "<<snp->GetLabel()<<" "<<score<<"\n";
		}
	}

	void Evaluate(const char *label, float newScore) {
		if (newScore > score ){ 
			score = newScore;
			if (this->snp)
				snp->ReduceInstanceCount();
			snp=NULL;
			label = label;
		}
	}			
	
	void Evaluate(SnpAligned *snp, float newScore) {
		if (newScore > score) {
			label = snp->GetLabel();
			score = newScore;	
			if (this->snp)
				snp->ReduceInstanceCount();
			this->snp=snp;
			if (snp)
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
	virtual void Append(float score, uint modelSize, uint test, const char *label)=0;
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
