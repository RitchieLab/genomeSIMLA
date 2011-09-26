//
// C++ Implementation: evalbalancedaccuracy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "evalbalancedaccuracy.h"

namespace ESE {

void EvalBalancedAccuracy::ReportResults(ostream *os) {
	*os<<"Quick cheesy test report:\nTop Model of A Balanced Accuracy Search:\n";
	uint count=results.size();

	for (uint i=0; i<count; i++) {
		if (distribution)
			*os<<"\t"<<results[i]->topModel->GetLabel()<<"\t"<<results[i]->highBA<<"\t"<<distribution->GetPValue(results[i]->topModel->GetLabelCount(), results[i]->highBA)<<"\n";
		else
			*os<<"\t"<<results[i]->topModel->GetLabel()<<"\t"<<results[i]->highBA<<"\n";
	}
	
}

}
