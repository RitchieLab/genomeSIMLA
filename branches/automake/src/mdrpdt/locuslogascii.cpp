//
// C++ Implementation: locuslog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locuslogascii.h"
#include "genetics/genetics.h"

namespace MDR {

using namespace Genetics;

void LocusLogAscii::WriteSnp(SnpAligned *snp) {
	*stream<<snp->GetLineNumber()<<" "<<snp->GetLabel()
		<<" "<<snp->GetGenotypeAlleles(0)
		<<" "<<snp->GetGenotypeAlleles(1)
		<<" "<<snp->GetGenotypeAlleles(2)
		<<" "<<snp->GetInclusionStatus()<<"\n";
}


}
