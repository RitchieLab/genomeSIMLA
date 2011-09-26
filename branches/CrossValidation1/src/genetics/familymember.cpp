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
	




}
