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
//	stringstream ss;
	

	uint localSize=genotypePresent[0].individuals.size();
	BitSetType missingData = other->genotypePresent[0].individuals;
/*	ss<<"  Size: "<<localSize<<"\n";
	ss<<"  Missing 1: "<<genotypePresent[0].individuals<<"\n";
	ss<<"  Missing 2: "<<other->genotypePresent[0].individuals<<"\n";
*/
	uint remoteSize=missingData.size();

	assert(localSize > 0);

	if (localSize < remoteSize) {
		for (uint each=0; each<genotypePresent.size(); each++) {
			genotypePresent[each].individuals.resize(remoteSize);
		}
	}
	else if (remoteSize < localSize) 
		missingData.resize(localSize);

	//missingData = ~missingData; 

	GenotypeData *gt=new GenotypeData(genoLookup);

	BitSetType dataPresent;
	dataPresent.resize(localSize);

	for (uint i=0; i< genotypePresent.size(); i++) {		
		dataPresent |= genotypePresent[i].individuals;
		BitSetType curGt = genotypePresent[i].individuals & missingData;
/*		ss<<"    "<<genotypePresent[i].individuals<<"\n";
		ss<<"    "<<missingData<<"\n";
		ss<<"    "<<curGt<<"\n\n";
*/
		gt->genotypePresent.push_back(GtStorage(genotypePresent[i].label.c_str(), curGt));
	}
	
	//cout<<"Masking two genotypes:\n"<<ss.str()<<"    "<<dataPresent<<"\n";
	assert(gt->genotypePresent.size() == genotypePresent.size());
	return gt;
}


int GenotypeData::GetGenotypeIndex(uint position) {
	if (genotypePresent.size() == 0)
		return 0;

	if (genotypePresent[0].individuals.size() <= position)
		return genoLookup->GetNotEncodedIdx();
	int count=genotypePresent.size();
	int idx=0;
	
	//cout<<"Looking for individual #"<<position<<" in the following data.\n";
	
	bool isFound = !genotypePresent[idx].individuals[position];

	for (idx=1; !isFound && idx<count; idx++) {
		//cout<<" "<<idx<<" "<<genotypePresent[idx].individuals[position]<<"\t"<<genotypePresent[idx].individuals<<"\n";
		isFound = isFound || genotypePresent[idx].individuals[position];
	}

	if (isFound)
		return idx -1;
	else
		return genoLookup->GetNotEncodedIdx();
}

void GenotypeData::ZeroGenotype( uint position) {
	VerifyDimensions(0, position+1);

	uint genotypeCount = genotypePresent.size();
	for (uint i=0; i< genotypeCount; i++) 
		genotypePresent[i].individuals[position]=i==0;
}


}
