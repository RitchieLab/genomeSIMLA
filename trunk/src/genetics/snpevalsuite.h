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
	/**
	 * @brief Main Constructor
	 * @param training This is the evaluation method method
	 * @param fDist The fitness distribution to be used
	 * @param orDist The odds ratio distribution to be used (if any)
	 */
	SnpEvalSuite(SnpEvaluation *training, PTestDistribution *fDist, PTestDistribution *orDist = NULL, PTestDistribution *peDist = NULL, PTestDistribution *lrDist = NULL);

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

	virtual void SetOddsRatioDist(PTestDistribution *dist);
	virtual void SetLRDist(PTestDistribution *dist);


	virtual PTestDistribution *GetOddsRatioDist();

	/**
	 * @brief Set the Distribution object to be used
	 */
	virtual void SetFitnessDist(PTestDistribution *dist);

	/**
	 * @brief Returns the distribution object used by the evaluation object
	 */
	virtual PTestDistribution *GetFitnessDist();

	virtual void SetPEDist(PTestDistribution *dist);
	virtual PTestDistribution *GetPEDist();
	virtual PTestDistribution *GetLRDist();

	void BuildDistribution();

	void PostEvaluation();

	SnpEvaluation *GetTrainer();
	ModelStatistics GetTopModel(uint modelSize, SnpAligned*& snp);

protected:
	SnpEvaluation *trainingSuite;
	vector<SnpEvaluation *> tests;
	PTestDistribution *fitnessDist;
	PTestDistribution *orDist;
	PTestDistribution *peDist;
	PTestDistribution *lrDist;

};

inline
ModelStatistics SnpEvalSuite::GetTopModel(uint modelSize, SnpAligned*& snp) {
	return trainingSuite->GetModelStats(modelSize, snp);
}

inline
void SnpEvalSuite::BuildDistribution() {
	uint count=tests.size();
	
	for (uint i=0; i<count; i++) {
		SnpAligned *topSnp;
		float averageTesting = tests[i]->GetAverageTesting(topSnp);
		fitnessDist->AppendTest(averageTesting, topSnp);
	}
	if (fitnessDist)
		fitnessDist->Sort();
}

inline
SnpEvaluation *SnpEvalSuite::GetTrainer() {
	return trainingSuite;
}

inline
PTestDistribution *SnpEvalSuite::GetFitnessDist() {
	return fitnessDist;
}

inline 
PTestDistribution *SnpEvalSuite::GetOddsRatioDist() {
	return orDist;
}

inline
PTestDistribution *SnpEvalSuite::GetPEDist(){ 
	return peDist;
}

inline
PTestDistribution *SnpEvalSuite::GetLRDist() {
	return lrDist;
}

inline
void SnpEvalSuite::SetPEDist( PTestDistribution *dist) {
	peDist = dist;
}

inline
void SnpEvalSuite::SetFitnessDist(PTestDistribution *dist) {
	fitnessDist=dist;
	trainingSuite->SetFitnessDist(dist);
}

inline
void SnpEvalSuite::SetOddsRatioDist(PTestDistribution *dist) {
	orDist = dist;
	trainingSuite->SetOddsRatioDist(dist);
}

inline
void SnpEvalSuite::SetTrainer(SnpEvaluation *method) {
	if (trainingSuite)
		delete trainingSuite;
	trainingSuite = method;
}

inline
void SnpEvalSuite::SetLRDist(PTestDistribution *dist) {
	lrDist = dist;
	trainingSuite->SetLRDist( dist);
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
SnpEvalSuite::SnpEvalSuite(SnpEvaluation *training, PTestDistribution *fDist, PTestDistribution *orDist, PTestDistribution *peDist, PTestDistribution *lrDist) : trainingSuite(training), fitnessDist(fDist), orDist(orDist), peDist(peDist), lrDist(lrDist) {
	training->SetFitnessDist( fDist );
	training->SetOddsRatioDist( orDist );
	training->SetPEDist( peDist );
}

inline
SnpEvalSuite::SnpEvalSuite() : trainingSuite(NULL), fitnessDist(NULL), orDist(NULL), peDist(NULL), lrDist(NULL) {
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
