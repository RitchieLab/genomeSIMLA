//
// C++ Implementation: minimalchromosome
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "minchrom.h"

namespace Simulation {

uint MinChrom::GetModelIndex(uint i) {
	assert(i < modelLoci->size());
	return (*modelLoci)[i];
}

uint MinChrom::GetModelSize() {
	return modelLoci->size();
}


}
