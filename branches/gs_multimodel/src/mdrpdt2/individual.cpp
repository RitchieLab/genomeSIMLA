//
// C++ Implementation: individual
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>
#include <iostream>
#include "utility/types.h"
#include "individual.h"

using namespace std;

namespace MdrPDT {


char Individual::AffectedValue		='2';
char Individual::UnaffectedValue	='1';

bool Individual::UseMerlin			=false;
Individual::Individual(const char *id, const char *pedID, int pedigreeIdx) : dropFromAnalysis(false), id(id), pedID(pedID), pat("0"), mat("0"), father(NULL), mother(NULL), pedigreeID(pedigreeIdx), status(-1), gender(0), validGenotypes(0) {
}


Individual::~Individual()	{
}


bool Individual::IsAffected() {
	return status;
}

void Individual::SetStatus(char status) {
	this->status = status == AffectedValue;

	//If status is unknown, then we shouldn't use this individual
	if (!this->status && status != UnaffectedValue) 
		DropFromAnalysis(true);
}

void Individual::SetGender(char gender) {
	this->gender=gender;
}

void Individual::IsAffected(bool isAffected) {
	status = isAffected;
}

void Individual::AddGenotype(char val) {
	assert(val!=0);
	genotypes.push_back(val);
	if (val > 0)
		validGenotypes++;
}

char Individual::GetVirtualGenotype(size_t idx) {
	assert(idx<virtualGenotypes.size());
	return virtualGenotypes[idx];
}

void Individual::SetGenotype(size_t idx, char gt) {
	assert(idx<genotypes.size());
	if (gt == 0 && genotypes[idx] > 0)
		validGenotypes--;
	genotypes[idx]=gt;
}
char Individual::GetGenotype(size_t idx) {
	assert(idx<genotypes.size());
	return genotypes[idx];
}
	
size_t Individual::CountGenotypes() {
	return genotypes.size();
}


void Individual::SetID(const char *id) {
	this->id = id;
}

string Individual::GetID() {
	return id;
}
void Individual::GetParentIDs(std::string& pat, std::string& mat) {
	pat = this->pat;
	mat = this->mat;
}
void Individual::SetParentIDs(const std::string& pat, const std::string& mat) {
	this->pat = pat;
	this->mat = mat;
}
bool Individual::HasVirtualGentoypes() {
	return virtualGenotypes.size() > 0;
}
int Individual::GetPedigreeID() {
	return pedigreeID;
}

void Individual::Write(ostream& os, bool reportVirtual) {
	string values[] = { "x x", "1 1", "1 2", "2 2" };
	std::vector<char >::iterator itr;
	std::vector<char >::iterator end;

	//Header stuff
	os<<pedID<<" "<<id<<" "<<pat<<" "<<mat<<" ";
	//6 vs 10 column header
	if (!UseMerlin) 
		os<<"0 0 0 "<<gender<<" 0 ";
	else
		os<<gender<<" ";
	

	//Since we don't want to analyze with it, but if it's to be treated as a parent, it's OK, we'll mark as unknown
	if (DropFromAnalysis()) 
		os<<"x";
	else if (IsAffected())
		os<<AffectedValue;
	else
		os<<UnaffectedValue;

	if (reportVirtual) {
		itr = virtualGenotypes.begin();
		end = virtualGenotypes.end();	
	} else {
		itr = genotypes.begin();
		end = genotypes.end();
	}
	while (itr != end) {
		int genotype=(int)*itr++;
		if (genotype>=0 && genotype<4)
			os<<" "<<values[genotype];
		else
			os<<" x x";
		}
}
void Individual::ReportGenotypes(ostream& os, bool reportVirtual, string sep) {
	std::vector<char >::iterator itr;
	std::vector<char >::iterator end;

	if (reportVirtual) {
		itr = virtualGenotypes.begin();
		end = virtualGenotypes.end();	
	} else {
		itr = genotypes.begin();
		end = genotypes.end();
	}

	while (itr != end) 
		os<<sep<<(int)*itr++;
	
}
void Individual::GenerateReport(ostream& os, bool reportVirtuals) {
	string status = "Unaffected ";
	if (IsAffected())
		status="Affected   ";
	os<<"\t"<<pedID<<"x"<<id<<"\t["<<pat<<":"<<mat<<"] "<<status<<": (";
	ReportGenotypes(os, false, "");
	
	os<<")\n";
	if (reportVirtuals) {
		os<<"           \t                       ~ : (";
		ReportGenotypes(os, true, "");
		os<<")\n";
	}

}

bool Individual::SetupVirtuals(ostream& os) {
	virtualGenotypes.clear();
	//For now, let's not mess with logging
	if (father == NULL || mother == NULL)
		return false;
	if (father->GetValidGenotypeCount() * mother->GetValidGenotypeCount() == 0) {
		os<<"Sibship: "<<id<<":"<<father->GetID()<<"x"<<mother->GetID()<<"\t\tUnable to create virtual sibling due to lack of parental genotypes\n";
		return false;
	}

	uint count = CountGenotypes();

	int f;
	int m;
	int c1;
	int badData=0;
	for (uint i=0; i<count; i++) {
		f=father->GetGenotype(i);
		m=mother->GetGenotype(i);
		c1=GetGenotype(i);

		//Let's make sure any parents that are missing data get filtered out as unknown
		if (f * m<=0  || c1<=0)
			virtualGenotypes.push_back(badData);
		//If the parents have the same genotype, we have two options. The children have same, unless both parents are heterozygote
		else if (m==f)
			if (m==2 && c1==1)
				virtualGenotypes.push_back(3);
			else if (m==2 && c1==3)
				virtualGenotypes.push_back(1);
			else if (c1 == f)					
				virtualGenotypes.push_back(c1);
			else {
//				os<<"-\tChild Data("<<i<<") Doesn't match parents: "<<c1<<" can't be produced from "<<m<<", "<<m<<"\n";
				virtualGenotypes.push_back(badData);
			}						
		else {
			//If the parents are seperated by 1, then we have a hetero/homo sitation. Technically, the children will mimic one or the other
			if (abs(m-f) < 2) {
				if (f == c1 || m == c1) {
					if (GetGenotype(i) == f)
						virtualGenotypes.push_back(m);
					else
						virtualGenotypes.push_back(f);
				}	
				else {
//					os<<"-\tChild Data("<<i<<") Doesn't match parents: "<<c1<<" can't be produced from "<<m<<", "<<f<<"\n";
					virtualGenotypes.push_back(badData);
				}						
			}
			//In this case, the parents have different homozygote genotypes. The children must be heterozygous
			else if (c1 == 2) 
				virtualGenotypes.push_back(2);
			else {
//				os<<"-\tChild data("<<i<<") doesn't match parents: "<<c1<<" can't be produced from "<<m<<", "<<f<<"\n";
				virtualGenotypes.push_back(badData);
			}
		}
	}
	return true;
}

}
