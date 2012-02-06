//
// C++ Implementation: casecontrolstatus
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "casecontrolstatus.h"

namespace Utility {


void StatusContainer::GenerateReport(ostream *os){
	size_t count=status.size();

	if (count > 0) {
		BitSetType totalAffected = status[0].affected;
		BitSetType totalUnaffected = status[0].unaffected;
		//Let's grab the first one and start there
		for (uint i=0; i<count; i++) {
			*os<<"Affected: "<<status[i].affected.count()<<"\tUnaffected: "<<status[i].unaffected.count()<<"\n";
			totalAffected = totalAffected | status[i].affected;
			totalUnaffected = totalUnaffected | status[i].unaffected;
		}
		*os<<totalAffected<<"\n"<<totalUnaffected<<"\n";	
	}
}	


}
