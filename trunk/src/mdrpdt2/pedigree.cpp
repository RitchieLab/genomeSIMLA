//
// C++ Implementation: pedigree
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigree.h"

using namespace std;

namespace MdrPDT {

Pedigree::Pedigree(const char *id, int pedID) : id(id), pedigreeID(pedID), dspCount(-1), dropFromAnalysis(false) {
}


Pedigree::~Pedigree()	{
	PurgeMembers();
}

void Pedigree::PurgeMembers() {
	std::map<std::string, Individual *>::iterator itr = members.begin();
	std::map<std::string, Individual *>::iterator end = members.end();
	
	while (itr != end) {
		delete itr->second;
		itr++;
	}
	
	members.clear();
}

int Pedigree::GetMemberCount() {
	return members.size();
}
Individual *Pedigree::GetMember(const char *id, bool createIfNotPresent) {
	Individual *ind = NULL;

	if (members.find(string(id)) == members.end()) {
		if (createIfNotPresent) {
			ind = new Individual(id, this->id.c_str(), pedigreeID);
			members[id] = ind;
		}
	} else  {
		ind = members[id];
	}
	return ind;
}

void Pedigree::InitSibships() {
	sibships.clear();
	//Work through the members, and create sibships
	//A sibship ONLY contains children. The parents will have to be obtained through the members[] structure
	std::map<std::string, Individual *>::iterator itr = members.begin();
	std::map<std::string, Individual *>::iterator end = members.end();
	if (DropFromAnalysis()) {
		return;
	}
	while (itr != end) {
		Individual *ind = itr->second;
		string mat, pat;
		ind->GetParentIDs(pat, mat);
		string parentIDs = pat + "x" + mat;

		if (parentIDs != "0x0" && !ind->DropFromAnalysis())
			sibships[parentIDs].AddIndividual(ind);
		itr++;	
	}
}


int Pedigree::WriteGenotypeData(char *data, int *meta, char *folds, int offset, Utility::Random& gen, int xvSlice, int stride, bool permute) {
	std::map<std::string, Sibship >::iterator itr = sibships.begin();
	std::map<std::string, Sibship >::iterator end = sibships.end();

	int members = offset;
	while (itr != end) {
		members=itr->second.WriteGenotypeData(data, meta, folds, members, gen, xvSlice, stride, permute);
		itr++;
	}
	return members;
}

int Pedigree::GetDSPCount() {
	if (dspCount <0) {
		std::map<std::string, Sibship >::iterator itr = sibships.begin();
		std::map<std::string, Sibship >::iterator end = sibships.end();
		dspCount = 0;
		while (itr != end) {
			dspCount += itr->second.CountDSPs();
			itr++;
		}
	}

	return dspCount;
}
void Pedigree::Reconcile() {
	std::map<std::string, Individual *>::iterator itr = members.begin();
	std::map<std::string, Individual *>::iterator end = members.end();

	while (itr!=end) {
		string fatherID, motherID;
		itr->second->GetParentIDs(fatherID, motherID);
		
		Individual *father = GetMember(fatherID.c_str(), false);
		Individual *mother = GetMember(motherID.c_str(), false);
		itr->second->SetFather(father);
		itr->second->SetMother(mother);

		//If we decide to keep up with the children associated with this parent, here is where we do it
		itr++;
	}
}
void Pedigree::ExpandSibships(ostream& os) {
	std::map<std::string, Sibship >::iterator itr = sibships.begin();
	std::map<std::string, Sibship >::iterator end = sibships.end();

	string fatherID, motherID;
	string sibshipID;
	//For now, we'll just expand trios
	while (itr != end) {
		Sibship &sibs = itr->second;
	
		if (sibs.GetParentIDs(fatherID, motherID)) {
			sibs.InitVirtualSibs(os);
		}
		int dspCount = sibs.GetDspCount();
		sibshipID = itr->first;
		itr++;
		if (dspCount < 1) 
			sibships.erase(sibshipID);

	}

	os<<"Pedigree: "<<id<<" Contributed "<<this->GetDSPCount()<<" DSPs.\n";

	
}

void Pedigree::Write(ostream& os) {
	std::map<std::string, Individual *>::iterator itr = members.begin();
	std::map<std::string, Individual *>::iterator end = members.end();
	if (DropFromAnalysis()) {
		return;
	}
	while (itr != end) {
		Individual *ind = itr->second;
		ind->Write(os);
		os<<"\n";
		itr++;
	}
}
void Pedigree::PostLoad(ostream& os) {
	InitSibships();

	//Work through each sibship, and expand trios (and possibly others as well)
	ExpandSibships(os);

	
}

}
