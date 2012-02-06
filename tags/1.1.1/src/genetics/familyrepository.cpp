//
// C++ Implementation: familyrepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "familyrepository.h"
#include "familyrepoevaluation.h"
#include <iomanip>
namespace Genetics {
using namespace std;

FamilyRepository *FamilyRepository::_instance = NULL;
int FamilyRepository::_instanceCount =0;

void FamilyRepository::PerformEvaluation( FamilyRepoEvaluation *eval) {
	map <string, FamilyNode *>::iterator itr = entries.begin();
	map <string, FamilyNode *>::iterator end = entries.end();

	FamilyNode *family;
	uint count=0;
	for (; itr != end; itr++) {
		family=itr->second;

		//Delete the whole family if the return tells you so
		if (eval->Evaluate(family, count++)) {
			//delete family;
			//entries.erase(itr);
		}
	}
}


void FamilyRepository::GenerateReport(ostream *os) {
	*os<<"Total Families Identified: "<<entries.size()<<"\n";
	*os<<"The following represents the contributions of each pedigree to the DSPs used in the analyses.\n";
	*os<<"Column 1 is the pedigree ID. 2 represents the total number of individuals contributed to the analysis. \n";
	*os<<setw(12)<<right<<"Pedidgree"<<setw(8)<<"Total"<<setw(8)<<"Unique"<<setw(8)<<"Unique"<<"   Sibship breakdown\n";
	*os<<setw(12)<<right<<"   ID    "<<setw(8)<<"Cont."<<setw(8)<<"Aff."<<setw(8)<<"Unaff"<<"      FxM [ sibs ] ...\n";
	if (os) {
		map <string, FamilyNode *>::iterator itr = entries.begin();
		map <string, FamilyNode *>::iterator end = entries.end();

		for (; itr!=end; itr++) 
			itr->second->GenerateReport(os);
	}
}	

}

