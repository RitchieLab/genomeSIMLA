//
// C++ Implementation: sibship
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "sibship.h"

using namespace std;

namespace MdrPDT {

Sibship::Sibship() : effectiveDspCount(-1)	{
}

Sibship::Sibship(const Sibship& other) {
	affectedMembers = other.affectedMembers;
	unaffectedMembers = other.unaffectedMembers;
	effectiveDspCount = other.effectiveDspCount;
}

Sibship::~Sibship(){
}


void Sibship::AddIndividual(Individual *ind) {
	if (ind->IsAffected())
		affectedMembers.push_back(ind);
	else
		unaffectedMembers.push_back(ind);

	//Reset the count, just in case it was evaluated already for some reason
	effectiveDspCount = -1;
}

int Sibship::CountDSPs() {
	effectiveDspCount = affectedMembers.size() * unaffectedMembers.size();
	
	vector<Individual*>::iterator itr = affectedMembers.begin();
	vector<Individual*>::iterator end = affectedMembers.end();
	
	while (itr != end) {
		Individual *ind = *itr;
		if (ind->HasVirtualGentoypes())
			effectiveDspCount+=1;
		itr++;	
	}
	return effectiveDspCount;
}

int Sibship::GetDspCount() {
	if (effectiveDspCount<0)
		return CountDSPs();
	else
		return effectiveDspCount;
}

void Sibship::Permute(vector<Individual*>& aff, vector<Individual*>& unaff, Utility::Random& gen) {
	vector<Individual*> allIndividuals = aff;
	allIndividuals.insert(allIndividuals.end(), unaff.begin(), unaff.end());

	int affCount = aff.size();
	int count = allIndividuals.size();
	aff.clear();
	unaff.clear();

	random_shuffle(allIndividuals.begin(), allIndividuals.end(), gen);
	
	for (int i=0; i<count; i++) {
		if (i<affCount)
			aff.push_back(allIndividuals[i]);
		else
			unaff.push_back(allIndividuals[i]);
	}
}

bool Sibship::GetParentIDs(string &father, string &mother) {
	bool isValid = false;

	if (affectedMembers.size() > 0) {
		affectedMembers[0]->GetParentIDs(father, mother);
		isValid = father != "0" && mother != "0";
	}
	return isValid;
}


void Sibship::InitVirtualSibs(ostream& os) {
	std::vector<Individual*>::iterator itr = affectedMembers.begin();
	std::vector<Individual*>::iterator end = affectedMembers.end();

	while (itr != end) {
		Individual *ind = *itr;
		ind->SetupVirtuals(os);
		itr++;
	}

	itr = unaffectedMembers.begin();
	end = unaffectedMembers.end();

	while (itr != end) {
		Individual *ind = *itr;
		ind->SetupVirtuals(os);
		itr++;
	}
}
char *Sibship::WriteRandomizedTrio(char *data, int *meta, char *folds, int offset, Individual *affected, int xvSlice, int stride) {
	char *pos = data+offset;			///<This is where we are, currently	
	stride--;							///<Don't want to overshoot (the one is the unaffected)
	int genoCount = affected->CountGenotypes();
	
	//PedigreeID needs to be an integer, which is is continuous, in a sense, so that there no (or very few) gaps
	//across all pedigrees. It's different from the real id...and is used to sum up the D-values correctly
	meta[offset] 	= affected->GetPedigreeID();
	meta[offset+1]  = affected->GetPedigreeID();
	
	folds[offset] 	= (char)xvSlice;
	folds[offset+1] = (char)xvSlice;

	for (int i=0; i<genoCount; i++) {		
		*(pos++) 		= affected->GetVirtualGenotype(i);		
		*(pos) 	= affected->GetGenotype(i);
		pos+=stride;
	}
	return pos++;
}
char *Sibship::WriteGenotypes(char *data, int *meta, char *folds, int offset, Individual *affected, int xvSlice, int stride) {
	//OK, we will save a little bit of time writing both the affected and it's mate out together
	char *pos = data+offset;			///<This is where we are, currently	
	stride--;							///<Don't want to overshoot (the one is the unaffected)

	int genoCount = affected->CountGenotypes();
	
	//PedigreeID needs to be an integer, which is is continuous, in a sense, so that there no (or very few) gaps
	//across all pedigrees. It's different from the real id...and is used to sum up the D-values correctly
	meta[offset] 	= affected->GetPedigreeID();
	meta[offset+1]  = affected->GetPedigreeID();
	
	folds[offset] 	= (char)xvSlice;
	folds[offset+1] = (char)xvSlice;

	for (int i=0; i<genoCount; i++) {
		*(pos++) = affected->GetGenotype(i);
		*(pos) = affected->GetVirtualGenotype(i);
		pos+=stride;
	}
	return pos++;
}

char *Sibship::WriteGenotypes(char *data, int *meta, char *folds, int offset, Individual *affected, Individual *unaffected, int xvSlice, int stride) {
	//OK, we will save a little bit of time writing both the affected and it's mate out together
	char *pos = data+offset;			///<This is where we are, currently
	stride--;							///<Don't want to overshoot (the one is the unaffected)

	int genoCount = affected->CountGenotypes();
	
	//PedigreeID needs to be an integer, which is is continuous, in a sense, so that there no (or very few) gaps
	//across all pedigrees. It's different from the real id...and is used to sum up the D-values correctly
	meta[offset] 	= affected->GetPedigreeID();
	meta[offset+1] 	= unaffected->GetPedigreeID();
	
	folds[offset] = (char)xvSlice;
	folds[offset+1] = (char)xvSlice;

	for (int i=0; i<genoCount; i++) {
		*(pos++) = affected->GetGenotype(i);
		*(pos) = unaffected->GetGenotype(i);
		pos+=stride;
	}
	return pos++;
}
int Sibship::WriteGenotypeData(char *data, int *meta, char *folds, int offset, Utility::Random& gen, int xvSlice, int stride, bool permute) {
	vector<Individual *> affecteds 		= affectedMembers;
	vector<Individual *> unaffecteds 	= unaffectedMembers;
	if (permute)
		Permute(affecteds, unaffecteds, gen);

	vector<Individual *> dsps;

	int affCount 	= affecteds.size();
	int unaffCount	= unaffecteds.size();
	int memberCount = CountDSPs() * 2;
	int members = offset;

	//If we are permuting, and the sibship has a single sibling, we have to handle permutation differently
	if (permute && affCount==1 && unaffCount==0) {
		if (gen.drand() < 0.5) 
			WriteRandomizedTrio(data,meta,folds,members,affecteds[0],xvSlice,stride);
		else
			WriteGenotypes(data, meta, folds, members, affecteds[0], xvSlice, stride);
		members+=2;
	} else {
		for (int a = 0; a<affCount; a++) {
			for (int u = 0; u<unaffCount; u++) {
				WriteGenotypes(data, meta, folds, members, affecteds[a], unaffecteds[u], xvSlice, stride);	
				members+=2;
			}
		}
		if (members < memberCount + offset) {
			//And do the non-transmitted DSPs
			for (int a=0; a<affCount; a++) {
				WriteGenotypes(data, meta, folds, members, affecteds[a], xvSlice, stride);
				members+=2;
			}
		}
	}
	return members;
}

}
