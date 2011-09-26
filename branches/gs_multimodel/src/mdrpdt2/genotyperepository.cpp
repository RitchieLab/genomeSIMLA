//
// C++ Implementation: genotyperepository
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genotyperepository.h"
#include <sstream>
#include <assert.h>

using namespace std;

namespace MdrPDT {

/**
 * Genotype Repository
 */
GenotypeRepository::GenotypeRepository() : genotypes(NULL), pedIDs(NULL), folds(NULL), rowWidth(0), indCount(0), snpCount(0) {

}
GenotypeRepository::GenotypeRepository(const GenotypeRepository& other) :
		rowWidth(other.rowWidth),
		indCount(other.indCount), 
		snpCount(other.snpCount)		 {
	cout<<"Ugh, why is this being copied?\n";
}
GenotypeRepository::~GenotypeRepository() {
	if (genotypes)
		delete[] genotypes;
	if (folds)
		delete[] folds;
	if (pedIDs) 
		delete[] pedIDs;
}	

//Eventually, we want might change these according to the locus
string GenotypeRepository::GetGenotypeEncoding(int snpIDx, int genotype) {
	static string gts[] = {"1/1", "1/2", "2/2"};
	assert(genotype<=2 && genotype>=0);
	return gts[genotype];
}
void GenotypeRepository::Dump(ostream& os) {
	cout<<"Individuals: "<<indCount<<"\n";
	cout<<"Snp Count  : "<<snpCount<<"\n";
	cout<<"             ---------------    The data\n";
	cout<<"Pedigree:\t";
	for (int i=0; i<indCount; i++) 
		cout<<" "<<pedIDs[i];
	cout<<"\n";
	cout<<"\nFolds:   \t";
	for (int i=0; i<indCount; i++) 
		cout<<" "<<(int)folds[i];
	cout<<"\n";
	cout<<"\nGenotypes:\n";
	for (int s=0; s<snpCount; s++) {
		cout<<"#"<<s<<"\t";
		for (int i=0; i<indCount; i++) {
			cout<<" "<<(int)genotypes[s*rowWidth+i];
		}
		cout<<"\n";
	}
}

void GenotypeRepository::InitExclusionList(vector<string>& list) {
	doAnalyzeSNP.resize(0);
	doAnalyzeSNP.resize(snpCount, true);
	vector<string>::iterator itr = list.begin();
	vector<string>::iterator end = list.end();

	
	while (itr != end){ 
		size_t idx = atoi(itr->c_str()) - 1;
		if (idx>=0 && idx<doAnalyzeSNP.size())
			doAnalyzeSNP[idx]=false;

		itr++;
	}
}


void GenotypeRepository::SetFoldCount(int foldCount) {
	dspCount.clear();
	dspCount.resize(foldCount, 0);
}
int GenotypeRepository::GetFoldCount() {
	return dspCount.size();
}
void GenotypeRepository::SetDSPCount(int fold,int count) {
	dspCount[fold]=count;
}
char *GenotypeRepository::Initialize(int snpCount, int indCount) {
	if (genotypes)
		delete[] genotypes;
	if (folds)
		delete[] folds;
	if (pedIDs) 
		delete[] pedIDs;

	//If we decide we want to add any extra columns, add them here
	rowWidth = indCount;
	this->indCount = indCount;
	this->snpCount = snpCount;
	
	genotypes 	= new char[snpCount * rowWidth];
	folds	  	= new char[rowWidth];
	pedIDs	 	= new int [rowWidth];
	//For now, I'm clearing the area. If we need some extra speed, we can undo this
	memset((void*)genotypes, 0, sizeof(char) * snpCount * rowWidth);
	memset((void*)folds, 0, sizeof(char) * rowWidth);
	memset((void*)pedIDs, 0, sizeof(int) * rowWidth);
	
	return genotypes;
}
int *GenotypeRepository::GetPedigreeIDs() {
	return pedIDs;
}
char *GenotypeRepository::GetFolds() {
	return folds;
}
char *GenotypeRepository::GetSNP(int idx) {
	assert(idx < snpCount);
	
	return genotypes +  (idx*rowWidth);
}

}
