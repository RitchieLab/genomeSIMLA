//
// C++ Implementation: snprepostxtsorted
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snprepostxtsorted.h"

namespace Genetics {
void SnpReposTxtSorted::Append(SnpAligned *newSnp) {
	ReportModel entry(newSnp);
	uint loci=newSnp->GetLabelCount();
	if (loci>contents.size())
		ResizeArray(loci);
	contents[loci-1]->push_back(entry);
	
}

/**
void SnpReposTxtSorted::Append(SnpAligned *newSnp) {
	char basicEntry[128];
	sprintf(basicEntry, "%s %d %s", newSnp->GetLineNumber().c_str(), newSnp->GetLastMdEval(), newSnp->GetTxtDescriptor().c_str());
	ReportEntry entry(newSnp->GetTxtDescriptor(), newSnp->GetLastMdEval(), basicEntry, basicEntry);
	uint loci=newSnp->GetLabelCount();
	if (loci>contents.size())
		ResizeArray(loci);
	contents[loci-1]->push_back(entry);
}
**/



}
