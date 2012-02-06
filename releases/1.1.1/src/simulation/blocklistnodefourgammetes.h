//
// C++ Interface: blocklistnodefourgammetes
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONBLOCKLISTNODEFOURGAMMETES_H
#define SIMULATIONBLOCKLISTNODEFOURGAMMETES_H
#include "blocklistnode.h"

namespace Simulation {


class BasicBlockFilter : public BlockFilter {
public:
	bool Evaluate(Locus &l) { return true; }

	BasicBlockFilter(): BlockFilter(true, 100) { }
	void Reset() { }
	bool IsValid() { return true; }
	int GetScore() { return score; }

};


/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BlockListNodeFourGammetes : public BlockListNode {
public:
	BlockListNodeFourGammetes(float criterion, BlockListNode *previous, BlockList *completedList, BlockList *deletedList);
    BlockListNodeFourGammetes(float criterion);
    virtual ~BlockListNodeFourGammetes();


	bool ValidateBlocks(vector<Locus>& loci, BlockFilter *filter);

	/**
	 * @brief generate summary report
	 */
	void Report(ostream& os);

	/**
	 * @brief Generate a more detailed report of the block contents
	 */
	void DetailedReport(ostream &os);

	virtual BlockListNode *Clone() {
		BlockListNode *copy = new BlockListNodeFourGammetes(criterion, previous, completedList, deletedList);
		((BlockListNodeFourGammetes *)copy)->avgMAF 		= avgMAF;
		((BlockListNodeFourGammetes *)copy)->blockDensity 	= blockDensity;
		((BlockListNodeFourGammetes *)copy)->firstIdx 		= firstIdx;
		((BlockListNodeFourGammetes *)copy)->lastIdx 		= lastIdx;
		((BlockListNodeFourGammetes *)copy)->firstLocation 	= firstLocation;
		((BlockListNodeFourGammetes *)copy)->lastLocation 	= lastLocation;
		((BlockListNodeFourGammetes *)copy)->loci 			= loci;
		((BlockListNodeFourGammetes *)copy)->maxAlleleFreq 	= maxAlleleFreq;
		((BlockListNodeFourGammetes *)copy)->minAlleleFreq	= minAlleleFreq;
		((BlockListNodeFourGammetes *)copy)->selectionScore	= selectionScore;
		((BlockListNodeFourGammetes *)copy)->totalDistance	= totalDistance;
		((BlockListNodeFourGammetes *)copy)->maxImageWidth 	= maxImageWidth;
		((BlockListNodeFourGammetes *)copy)->rowsWritten	= rowsWritten;
	return copy;
}
	
	/**
	 * @brief This is used for sorting in STD sort functions
	 * @note We want it to work out so the best block is the first in the vector
	 */
	bool operator<(const BlockListNodeFourGammetes& other) const;
	BlockListNode *NewNode(BlockListNode *previous, uint firstIdx, uint lastIdx);


	/**
	 * @Brief starts the summary (in HTML format)
	 */
	void HtmlSummaryBegin(ostream &os);	

	/**
	 * @brief write block details to the stream (HTML format)
	 * @param os Stream to be written to
	 * @param linkToDetails The name of the detailed plot to be imbedded in the report
	 */
	void HtmlSummary(ostream &os, const char *linkToDetails);
		
	/**
	 * @brief Close out the block's report
	 */
	void HtmlSummaryEnd(ostream &os);

	/**
	 * @brief Sets up the detailed hmtl report
	 * @return the name of the anchor to which the page can be linked
	 */
	string HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails);
protected:
	size_t rowsWritten;
	size_t maxImageWidth;
};

}

#endif
