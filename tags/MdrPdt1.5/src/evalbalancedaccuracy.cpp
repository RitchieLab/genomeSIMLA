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

namespace MDR {

void EvalBalancedAccuracy::ReportResults(ostream *os) {
	*os<<"Quick cheesy test report:\nTop Model of A Balanced Accuracy Search:\n";
	uint count=results.size();

	for (uint i=0; i<count; i++) {
		*os<<"\t"<<results[i]->topModel->GetLabel()<<"\t"<<results[i]->statistics.trainingFitness;
/**--This is for when we redo cross validation
		if (foldCount > 0)
			*os<<"\t"<<results[i]->statistics.testingFitness;
*/		if (fitnessDist)
			*os<<"\t"<<fitnessDist->GetPValue(results[i]->statistics.trainingFitness, results[i]->topModel->GetLabelCount())<<"\n";
	}
	
}

}
