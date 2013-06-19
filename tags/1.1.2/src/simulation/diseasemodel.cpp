//
// C++ Implementation: diseasemodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "diseasemodel.h"
namespace Simulation {
namespace StatusModel {

void DiseaseModel::AddDiseaseLoci(const char *label, uint chromID, uint locus, float al1) {
	loci.push_back(DiseaseLocus(label, chromID, locus, al1, 1.0-al1));
}


bool DiseaseModel::GetDiseaseLoci(int chromID, vector<uint> &relLoc) {
	relLoc.clear();
	
	bool anyFound = false;
	ModelLociArray::iterator itr = loci.begin();
	ModelLociArray::iterator end = loci.end();

	while (itr != end) {
		if (itr->chromosome == chromID) {
			anyFound = true;
			relLoc.push_back(itr->locusIdx);
		}
		itr++;
	}
	return anyFound;
}



DiseaseLocus DiseaseModel::GetLocus(uint pos) {
	return loci[pos];
}


BitSetType DiseaseModel::AssociatedChromosomes(uint chromCount) {
	BitSetType chrom(chromCount, false);
	if (chromCount) {
		ModelLociArray::iterator itr = loci.begin();
		ModelLociArray::iterator end = loci.end();
	
		while (itr != end) { 
			chrom[itr->chromosome] = true;
			itr++;
		}
	}
	return chrom;
}

}
}
