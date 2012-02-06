//
// C++ Implementation: locusreport
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locreport.h"
#include "chrompool.h"

namespace Simulation {

namespace Visualization {

/**
 * @note blocksOfConcern is a copy! I don't want to be sorting the original block list
 */
void LocusReport::Init(vector<Locus>& loci, BlockList *blocksOfConcern, const char *blockReport, const char *chromLabel) {
	vector<Locus>::iterator locItr = loci.begin();
	vector<Locus>::iterator locEnd = loci.end();

	BlockList::iterator blkItr = blocksOfConcern->begin();
	BlockList::iterator blkEnd = blocksOfConcern->end();

	//We have to sort this list so that it is in order by location
	sort(blkItr, blkEnd, SortBlocksByLocation());

	while (locItr != locEnd) {
		vector<LocusSelection>::iterator selItr = selectors.begin();
		vector<LocusSelection>::iterator selEnd = selectors.end();

		string locID = locItr->GetLabel();
		while (selItr != selEnd) {
			bool continueToSearchBlocks = selItr->EvaluateRange(*locItr);

			if (continueToSearchBlocks) {
				LocusReportEntry *entry = new LocusReportEntry();
				entry->locus = *locItr;
				entry->blockDescriptionFile = blockReport;
				entry->chromosomeID  = chromLabel;
				uint location = locItr->GetLocation();

				blkItr = blocksOfConcern->begin();
				blkEnd = blocksOfConcern->end();

				while (blkItr != blkEnd && continueToSearchBlocks) {
					//We are working on the assumption that blocks are sorted by the first SNP's position. 
					//Thus, once we find a place where the SNP is too far to the left, we can move on
					if (location >= (*blkItr)->GetFirstLocation()) {
						if (location <= (*blkItr)->GetLastLocation()) {
							entry->associatedBlocks.push_back((*blkItr)->Clone());
						}
					}
					else
						continueToSearchBlocks = false;
					
					blkItr++;
				}
				if (! (selItr->Evaluate(entry) > 0.0000))
					delete entry;
			}

			selItr++;
		}
		locItr++;
	}

}

void LocusReport::Reset()	{
	vector<LocusSelection>::iterator itr = selectors.begin();
	vector<LocusSelection>::iterator end = selectors.end();

	while (itr != end) {
		itr->Reset();
		itr++;
	}
}


void LocusReport::WriteHtmlReport(ostream& report) {
	report<<"<LINK REL='stylesheet' TYPE='text/css' HREF=\""<<Simulation::ChromPool::cssFilename<<"\">\n<DIV class='image-frame'>\n";

	vector<LocusSelection>::iterator itr = selectors.begin();
	vector<LocusSelection>::iterator end = selectors.end();

	while (itr != end) {
		itr->WriteHtmlReport(report);
		itr++;
	}

}

LocusReport::~LocusReport() {
}


}

}
