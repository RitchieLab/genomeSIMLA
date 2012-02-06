//
// C++ Interface: snpevalsuite
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONSNPEVALSUITE_H
#define GENETICS_EVALUATIONSNPEVALSUITE_H

#include "snpevaluation.h"

namespace Genetics {

namespace Evaluation {


/**
@brief Provides structure for evaluating SNP models (including x-validation and PTests)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpEvalSuite : public Genetics::ValidationFunctors::SnpVerificationMethod {
public:
    SnpEvalSuite();
	SnpEvalSuite(SnpEvaluation *training, PermutationTestDist *dist);

    ~SnpEvalSuite();

	/**
	 * @brief Set the training set 
	 */
	void SetTrainer(SnpEvaluation *);

	/**
	 * @brief Add a set to the permutation vector
	 */
	void AddPTest(SnpEvaluation *);
	/**
	 * @brief Evaluate a model over the ptests and the status set. 
	 * This represents a single evaluation as opposed to a x-validation one	
	 */
	bool EvaluateModel(SnpAligned *snp, uint position);

	/**
	 * @brief Write the results to the desired stream
	 */  
	void ReportResults(ostream *os);
	void ReportPTests(ostream *os);


	/**
	 * @brief Set the Distribution object to be used
	 */
	virtual void SetDistribution(PermutationTestDist *dist);

	/**
	 * @brief Returns the distribution object used by the evaluation object
	 */
	virtual PermutationTestDist *GetDistribution();

	void BuildDistribution();

	void PostEvaluation();

protected:
	SnpEvaluation *trainingSuite;
	vector<SnpEvaluation *> tests;
	PermutationTestDist *distribution;

};

inline
void SnpEvalSuite::BuildDistribution() {
	uint count=tests.size();
	
	for (uint i=0; i<count; i++) {
		SnpAligned *topSnp;
		float topScore = tests[i]->GetTopScore(topSnp);
		distribution->Append(topScore, i, topSnp);
	}
	if (distribution)
		distribution->Sort();
}

inline
PermutationTestDist *SnpEvalSuite::GetDistribution() {
	return distribution;
}

inline
void SnpEvalSuite::SetDistribution(PermutationTestDist *dist) {
	distribution=dist;
	trainingSuite->SetDistribution(dist);
}

inline
void SnpEvalSuite::SetTrainer(SnpEvaluation *method) {
	if (trainingSuite)
		delete trainingSuite;
	trainingSuite = method;
}

inline
void SnpEvalSuite::AddPTest(SnpEvaluation *test) {
	tests.push_back(test);
}


inline
void SnpEvalSuite::ReportResults( ostream *os) {
	//This will probably require some preparation type activities 
	// * Create distribution of PTests
	// * Prep report with p-test Results
	trainingSuite -> ReportResults( os );
}

inline
void SnpEvalSuite::ReportPTests( ostream *os ){
	//This will probably act on the distribution that was created
	cout<<"PTest report has not be configured.\n";
	assert(0);
}

inline
bool SnpEvalSuite::EvaluateModel(SnpAligned *snp, uint position) {
	snpsSeen++;
	bool success = trainingSuite->EvaluateModel(snp, position);

	uint count=tests.size();
	for (uint i=0; i<count; i++) 
		tests[i]->EvaluateModel(snp, position);	

	return success;
}


inline
void SnpEvalSuite::PostEvaluation() {
	trainingSuite->PostEvaluation();
}


inline
SnpEvalSuite::SnpEvalSuite(SnpEvaluation *training, PermutationTestDist *dist) : trainingSuite(training), distribution(dist) {
	training->SetDistribution( dist );
}

inline
SnpEvalSuite::SnpEvalSuite() : trainingSuite(NULL), distribution(NULL) {
}

inline
SnpEvalSuite::~SnpEvalSuite() {
	if (trainingSuite)
		delete trainingSuite;
	
	uint count=tests.size();
	for (uint i=0; i<count; i++) 
		delete tests[i];
	tests.clear();

}




}

}

#endif
