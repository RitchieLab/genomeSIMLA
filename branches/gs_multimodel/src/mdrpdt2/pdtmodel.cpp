/**
 * PdtModel 
 * 
 * Description:		This class encapsulates basic model details required for evaluating a model
 * 
 * Author: 		Eric Torstenson, (C) Marylyn Ritchie 2008
 * 
 * Copyright:		See COPYING file that comes with this distribution.
 */
#include "pdtmodel.h"
#include <vector>

using namespace std;
namespace MdrPDT {

int PdtModel::MaxLabelLength = 0;
using namespace Evaluation;

PdtModel::PdtModel() : modelID(""), foldData(NULL), pedigreeData(NULL), xvCount(0), pedCount(0), individualCount(0){ }

PdtModel::~PdtModel() { }

bool PdtModel::operator<(const PdtModel& other) { return GetTrainingT(0) < other.GetTrainingT(0); } 

void PdtModel::ReportFold(int fold, std::ostream& os) {
	testing[fold].Report(os);
}

float PdtModel::GetTrainingT(int fold) const {
	return training[fold].GetTStatistic();
}

float PdtModel::GetTestingT(int fold) const {
	return testing[fold].GetTStatistic();
}

float PdtModel::EvaluateFold(float &training, float &testing, int fold) {
	training = this->training[fold].GetTStatistic();
	testing = this->testing[fold].GetTStatistic();
	return testing;
}

void PdtModel::AddSnp(int snpID, char *ptr) {
	stringstream ss;
	ss<<modelID;
	if (MaxLabelLength > 0) 
		ss<<setw(MaxLabelLength)<<(*snpLabels)[snpID-1]<<" ";
	else if (modelID.length() > 0) 
		ss<<" "<<(*snpLabels)[snpID-1];
	else 
		ss<<(*snpLabels)[snpID-1];
	
	snpIDs.push_back(snpID);
	modelID = ss.str();
	loci.push_back(ptr);
}

void PdtModel::AddFold(tcalc& training, tcalc& testing) {
	this->training.push_back(training);
	this->testing.push_back(testing);
}

std::string PdtModel::GetPrintableLabel() {
	return "[ "+modelID+" ]";
}

std::string PdtModel::GetModelID() {
	return modelID;
}

int PdtModel::ModelSize() { 
	return loci.size();
}

int PdtModel::GetLocusID(int idx) {
	return snpIDs[idx];
}

char *PdtModel::GetSnpData(int idx) {
	return loci[idx];

}

void PdtModel::InitModelCounts() {
	//First, we need to gather the details we are about to use
	int modelSize			= snpIDs.size();
	int foldCount			= xvCount;
	vector<char *> genotypes = loci;
	
	makeGenotype = GenotypeConverter(modelSize, individualCount);
	genotypeCounter = OverallCounter(foldCount, pedCount, makeGenotype.genotypeCount);
	int *pedID				= pedigreeData;
	char *folds				= foldData;
	
	int pedigree, fold;
	int aff=0, unaff=0;
	int indCount = 0;
	
	//Figure out who has which counts
	while (indCount < individualCount) {
		pedigree = (uint)*pedID;
		fold = (uint)*folds;
		
		//Make sure we have valid genotypes
		if (makeGenotype.GetGenotype(genotypes, indCount, aff, unaff)) 
			genotypeCounter.AddCounts(fold, pedigree, aff, unaff);

		//Move forward by two
		pedID+=2;
		folds+=2;
		indCount+=2;
	}
}
MatchedOddsRatio PdtModel::EvaluateMOR() {
	return genotypeCounter.EvaluateMOR();
}
MatchedOddsRatio PdtModel::EvaluateMOR(int fold) {
	return genotypeCounter.EvaluateMOR(fold);
}
MatchedOddsRatio PdtModel::EvaluateBiasedMOR(int fold) {
	return genotypeCounter.EvaluateBiasedMOR(fold);
}
void PdtModel::InitPedigree(int indCount, int count, int *data) { 
	pedigreeData = data;
	pedCount = count;
	individualCount = indCount;
}

void PdtModel::CrossValidations(int count, char *data) {
	foldData = data;
	xvCount  = count;
}

char *PdtModel::FoldData()  {
	return foldData;
}

int *PdtModel::PedigreeData() {
	return pedigreeData;
}

int PdtModel::FoldCount() { 
	return xvCount;
}

int PdtModel::PedigreeCount() {
	return pedCount;
}

int PdtModel::IndividualCount() {
	return individualCount;
}

vector<int> PdtModel::GetModelIDs() { 
	return snpIDs;
}

int PdtModel::GetModelSize() {
	return snpIDs.size();
}

string PdtModel::GetSnpLabel(int n) {
	return (*snpLabels)[snpIDs[n]-1];
}
int PdtModel::GetSnpID(int n) {
	return snpIDs[n];
}
}
