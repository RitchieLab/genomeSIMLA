//
// C++ Implementation: snprepostxtsorted
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "powerrepos.h"

namespace MDR {

namespace Power {

void PowerRepos::Append(SnpAligned *newSnp) {
	ReportModel entry(newSnp);
	contents.push_back(entry);
	
}





}


}
