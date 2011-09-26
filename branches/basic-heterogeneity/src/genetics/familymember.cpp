//
// C++ Implementation: familymember
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "familymember.h"

namespace Genetics {
FamilyMemberPool *FamilyMemberPool::_instance = NULL;
int FamilyMemberPool::_instanceCount = 0;

FamilyMember::GenoLkup FamilyMember::lkupTable;
int FamilyMember::_affectedValue 		= 2;
int FamilyMember::_unaffectedValue 		= 1;
	


ostream *FamilyMember::Report( ostream *os) {
	string m = "0";
	string d = "0";
	
	if (mother)
		m=mother->GetID();
	if (father)
		d=father->GetID();
	string status = "Unaffected ";
	if (IsAffected())
		status="Affected   ";
	else if (IsUnknownStatus())
		status="Unknown   ";
	*os<<"\t"<<famID<<"x"<<id<<"\t["<<d<<":"<<m<<"] "<<status<<": (";
	if (genotypedata)
		genotypedata->Report(os);
	else 
		*os<<"no genotype data available";
	*os<<")\n";

	return os;
}


}
