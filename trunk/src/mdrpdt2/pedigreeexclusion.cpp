//
// C++ Implementation: pedigreeexclusion
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreeexclusion.h"
#include <sstream>
#include <iomanip>


using namespace std;

namespace MdrPDT {

namespace Validator {

PedigreeExclusion::PedigreeExclusion(vector<string>& pedigrees) : 
		exclusions(pedigrees.begin(), pedigrees.end()), exclusionCount(0), excludedPedigrees("")	{
}


PedigreeExclusion::~PedigreeExclusion()
{
}

void PedigreeExclusion::Evaluate(Pedigree* ped) {
	if (find(exclusions.begin(), exclusions.end(), ped->ID()) != exclusions.end()) {
		ped->DropFromAnalysis();
		excludedPedigrees+=" " + ped->ID();
		exclusionCount++;
	}
}

}

}
