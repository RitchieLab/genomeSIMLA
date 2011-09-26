//
// C++ Implementation: snp
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genotypedata.h"

namespace Genetics {

GenotypeData *GenotypeData::BuildMaskedGT(GenotypeData *other) {
	uint localSize=genotypePresent[0].individuals.size();
	BitSetType missingData = other->genotypePresent[0].individuals;
	uint remoteSize=missingData.size();

	assert(localSize > 0);

	if (localSize < remoteSize) {
		for (uint each=0; each<genotypePresent.size(); each++) {
			genotypePresent[each].individuals.resize(remoteSize);
		}
	}
	else if (remoteSize < localSize) 
		missingData.resize(localSize);

	missingData = ~missingData; 

	GenotypeData *gt=new GenotypeData(genoLookup);

	for (uint i=0; i< genotypePresent.size(); i++) {		
		BitSetType curGt = genotypePresent[i].individuals & missingData;
		gt->genotypePresent.push_back(GtStorage(genotypePresent[i].label.c_str(), curGt));
	}
	assert(gt->genotypePresent.size() == genotypePresent.size());
	return gt;
}


int GenotypeData::GetGenotypeIndex(uint position) {
	if (genotypePresent.size() == 0)
		return 0;

	if (genotypePresent[0].individuals.size() <= position)
		return genoLookup->GetNotEncodedIdx();
	bool isFound=false;
	int count=genotypePresent.size();
	int idx;
	for (idx=0; !isFound && idx<count; idx++) 
		isFound = isFound || genotypePresent[idx].individuals[position];

	if (isFound)
		return idx -1;
	else
		return genoLookup->GetNotEncodedIdx();
}

}
