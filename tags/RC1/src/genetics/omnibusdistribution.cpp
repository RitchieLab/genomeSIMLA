//
// C++ Implementation: omnibusdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "omnibusdistribution.h"

namespace Genetics {

namespace Reporting {



//Assumes that the vector has been sorted
//Notice that the modelsize is passed, even though omnibus doesn't care about that. That is just
//to satisfy the inheritance
float OmnibusDistribution::GetPValue(uint modelsize, float score) {
	assert(isSorted); 

	uint count=scores.size();
	float pValue=1.0;

	for (uint i=count; i>0; i--){ 
		if (scores[i-1].score<score) {		
			pValue=(float)(count-i+1)/(float)count;
			break;
		}
	}
	return pValue;
}


}

}
