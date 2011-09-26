//
// C++ Interface: ldplotter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLDPLOTMGR_H
#define SIMULATIONLDPLOTMGR_H
//#include <pngwriter.h>
#include <vector>

#include "locus.h"
#include "chromosome.h"
#include "ldwriter.h"
#include "blocklistnode.h"

#include "blocklistnodefourgammetes.h"
#include "simulation.h"
#include "ldcalculator.h"

namespace Simulation {

namespace Visualization {

using namespace std;


/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
	@brief Perform LD Analysis and write the plots out according to user's configuration
	
	
*/
class LdPlotManager{
public:

	/**
	 * @brief Construction
	 * @param results This is the ld data that will be written to file
	 */
    LdPlotManager(vector<Locus> * loci);

    ~LdPlotManager();

	/**
	 * @brief Setup the LD values.
	 * 
	 */
	void Init(vector<vector<LdResult> >& results, BlockListHead<BlockListNodeFourGammetes>* haplotypeBlocks);
	

	/**
	 * @brief Produce a DPrime ld plot
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @param imgHeight Returns the actual height of the chart
	 * @param imgWidth Returns the actual width of the chart
	 * @param ldBufferSize The number of snps on either side of the block that will show up in the plot
	 * @param Lbl Text label used as part of the final filename
	 * @return Actual filename produced
	 */
	string WriteLdDPrime(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl);
	//string WriteLdDPrime(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl);

	/**
	 * @brief Produce a ld plot based on rsquared
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @param imgHeight Returns the actual height of the chart
	 * @param imgWidth Returns the actual width of the chart
	 * @param ldBufferSize The number of snps on either side of the block that will show up in the plot
	 * @param Lbl Text label used as part of the final filename
	 * @return Actual filename produced
	 */
	string WriteLdRSquared(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl);
	//string WriteLdRSquared(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl);

	/**
	 * @brief Write the LD values to a plain text file	 
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @return Actual filename produced
	 */
	string WriteLdReport(const char *filename, BlockListNode *block);
	//string WriteLdReport(const char *filename, HaplotypeBlock *block);

	/**
	 * @brief Generic LD production routine- ld values for region surrounding the block are calculated (or reused)
	 * @param filename Basis for the actual filename to be used
	 * @param block The block around which the plot will be built
	 * @param writer The object that will actually render the desired details (text or graphical)
	 * @param padding The number of SNPs on either side of a selected block (if a block is selected). Otherwise, this has no effect
	 * @return Actual filename used
	 */
	string WriteLdData(const char *filename, BlockListNode *block, LdWriter* writer, uint padding);
	//string WriteLdData(const char *filename, HaplotypeBlock *block, LdWriter* writer, uint maxSnpDepth);

	/**
	 * @brief Generic LD production routine- ld values for region surrounding the block are calculated (or reused)
	 * @param filename Basis for the actual filename to be used
	 * @param first The index of the first snp to be rendered
	 * @param last The index of the last snp to be rendered
	 * @param writer The object that will actually render the desired details (text or graphical)
	 * @param maxSnpDepth The most snps a single snp can be paired with
	 * @return Actual filename used
	 */
	string WriteLdData(const char *filename, uint first, uint last, LdWriter* writer, uint maxSnpDepth);

	static uint maxSnpDistance;			///<The max number of basepairs allowed for LD Consideration

	
	/**
	 * @brief Set a label for the chart
	 */
	void SetLabel(const char *lbl);
	

	void CountRC(size_t first, size_t last, size_t &rows, size_t &cols, int &firstLocation, int &lastLocation);
protected:		
	/**
	 * @brief Vector of haplotype blocks that will be used by client objects to perform certain tasks
	 */
	BlockListHead<BlockListNodeFourGammetes> *haplotypes;
	BlockList *blocks;
	

	size_t minBlockDensity;				///<Min. Snp density seen for the region, used for calculating the SNP density bars
	size_t maxBlockDensity;				///<Max. Snp density seen for the region, used for calculating the SNP density bars
	string label;						///<Label that goes onto the plot
//	vector<LdResult>* results;			///<Array of statistics
	size_t *ldStart;					///<The index into the array for a given snp (this array is 1 entry 
										///<larger than the number of legal sub arrays

	vector<vector<LdResult> >* results;
	vector<Locus> * loci;
};
}
}

#endif
