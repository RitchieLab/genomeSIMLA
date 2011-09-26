//
// C++ Implementation: familyrepository
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "familyrepository.h"
#include <sstream>
#include <iomanip>
#include "foldproduction.h"
#include "pdtmodel.h"

using namespace std;

namespace MdrPDT {

bool PedigreeRepository::verboseFoldingReport = false;
/**
 * Family Repository
 */
PedigreeRepository::PedigreeRepository(): dspCount(0), snpCount(-1) { 
	
}

PedigreeRepository::~PedigreeRepository() {
	std::map<std::string, Pedigree*>::iterator itr = pedigrees.begin();
	std::map<std::string, Pedigree*>::iterator end = pedigrees.end();

	while (itr != end) {
		delete itr->second;
		itr++;
	}
}

void PedigreeRepository::LoadDat(const char *filename) {
	ifstream file(filename);
	snpLabels.clear();
	uint widestLabel = PdtModel::MaxLabelLength;
	while (!file.eof()) {
		string descriptor="", label="";
		file>>descriptor>>label;

		if (label.length() > widestLabel)
			widestLabel = label.length();
		
		if (label.length()>0 && descriptor=="M")
			snpLabels.push_back(label);
	}
	PdtModel::MaxLabelLength = widestLabel;
	assert(snpCount == (int)snpLabels.size());
}

bool PedigreeRepository::InitializeData(GenotypeRepository& genotypeData, Utility::Random& rnd, int xvCount, bool randomizeStatus) {
	char *data = genotypeData.Initialize(snpCount, dspCount * 2);
	int  *meta = genotypeData.GetPedigreeIDs();
	char *folds = genotypeData.GetFolds();

	//For now, I'm letting the families write the data....it's seems a little fishy, though
	std::map<std::string, Pedigree*>::iterator itr = pedigrees.begin();
	std::map<std::string, Pedigree*>::iterator end = pedigrees.end();
		
	//If we actually start using any extra stuff at the front of a row, we need to account for that
	int stride = dspCount * 2;

	//We are assuming 1 cross validation is here. I am guessing we need 
	//to make that a special function call
	int xv = 0;
	int *pedID = meta;
	int members = 0;
	while (itr != end) {
		members= itr->second->WriteGenotypeData(data, pedID, folds, members, rnd, xv, stride, randomizeStatus);
		itr++;
	}

	//Perform Fold generation
	FoldProduction foldProduction(xvCount);
	vector<FamilyRecord> records;
	int totalDSPs = BuildFamilyRecords(records);

	foldProduction.BuildXVFolds(folds, records, totalDSPs, rnd);
	if (!randomizeStatus && PedigreeRepository::verboseFoldingReport)
		foldProduction.GenerateReport(cout);
	
	genotypeData.SetFoldCount(xvCount);
	for (int i=0; i<xvCount; i++) {
		int count = foldProduction.GetDSPCount(i);
		genotypeData.SetDSPCount(i, count);
	}
	genotypeData.SetSnpLabels(&snpLabels);
	return true;
}

int PedigreeRepository::BuildFamilyRecords(vector<FamilyRecord>& records) {
	std::map<std::string, Pedigree*>::iterator itr = pedigrees.begin();
	std::map<std::string, Pedigree*>::iterator end = pedigrees.end();

	records.clear();
	int offset = 0;
	while (itr != end) {
		Pedigree *family = (itr++)->second;
		if (!family->DropFromAnalysis()) {
			int localOffset = family->GetDSPCount() * 2;
			records.push_back(FamilyRecord(offset, localOffset, family->ID().c_str()));
			offset+=localOffset;
		}
	}
	return offset;
}


uint PedigreeRepository::PostLoad(ostream& os) {
	//Finish up the details
	std::map<std::string, Pedigree*>::iterator itr = pedigrees.begin();
	std::map<std::string, Pedigree*>::iterator end = pedigrees.end();
	os<<"\n\n-------------------- DSP Contribution ------------------";
	dspCount = 0;
	while (itr != end) {
		itr->second->PostLoad(os);
		dspCount += itr->second->GetDSPCount();
//		itr->second->PostLoad();
		itr++;
	}
	if (dspCount < 1) {
		cerr<<"There are not enough affected individuals present to perform an analysis. Please check that affected status is designated according to: A)"
			<<Individual::AffectedValue<<"\tU)"<<Individual::UnaffectedValue<<". See the manual for instructions on how to correctly designate affected status using the configuration file.\n";
		exit(1);
	}
	return dspCount;
}
uint PedigreeRepository::Load(const char *filename, ostream& os) {
	Utility::LineParser lp('#');
	uint lines = lp.Parse(filename, this);	
	Iterator itr = GetIterator();
	Pedigree *ped = itr.GetNext();

	while (ped) {
		ped->Reconcile();
		ped = itr.GetNext();
	}

	char label[64];
	uint maxWidth = PdtModel::MaxLabelLength;
	//Unless we have a dat file, we'll just keep a simple integer based id
	for (int i=0; i<snpCount; i++) {
		sprintf(label, "%d", i+1);
		snpLabels.push_back(label);
		if (strlen(label) > maxWidth)
			maxWidth = strlen(label);
	}
	PdtModel::MaxLabelLength = maxWidth;
	return lines;
}

void PedigreeRepository::Write(const char *filename) {
	Iterator itr=GetIterator();
	ofstream file(filename);
	Pedigree *family = itr.GetNext();
	while (family) {
		if (!family->DropFromAnalysis()){ 
			family->Write(file);
		}
		family=itr.GetNext();
	}
}

Pedigree *PedigreeRepository::GetPedigree(const char *pedID, bool createIfNotPresent) {
	Pedigree *ped = NULL;
	
	if (pedigrees.find(pedID) != pedigrees.end())  {
		ped = pedigrees[pedID];
	}
	else if (createIfNotPresent) {
		ped = new Pedigree(pedID, pedigrees.size());
		
		pedigrees[pedID] = ped;
	}
	return ped;
}


Individual *PedigreeRepository::GetIndividual(const char *pedID, const char *indID, bool createIfNotPresent) {
	Individual *ind = NULL;
	Pedigree *family = NULL;
	
	family = GetPedigree(pedID, createIfNotPresent);
	if (family) 
		ind = family->GetMember(indID, createIfNotPresent);	
	
	return ind;
}



}
