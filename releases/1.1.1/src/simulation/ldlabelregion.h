//
// C++ Interface: ldlabelregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDLABELREGION_H
#define SIMULATION_VISUALIZATIONLDLABELREGION_H

#include "ldpngcomponent.h"
#include "locus.h"
//#include "haplotypeblock.h"
#include "blocklistnodefourgammetes.h"

namespace Simulation {

namespace Visualization {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

class LdLabelRegion : public LdPngComponent {
public:
    LdLabelRegion(ImageParameters *param);
    virtual ~LdLabelRegion();

	/**
	 * @brief Renders a locus' label to the header portion of the graph, if possible
	 */
	void AddLocus(Locus &l);

	/**
	 * @brief Perform any initial rendering (like the guides)
	 */
	void Open(const char *name);

	/**
	 * @brief Do any last minute cleanup or rendering
	 */
	void Close();

	/**
	 * @brief Determines a handful of properties to be used later
	 * @param yOffset indicates where the top of the label section is
	 * @param firstIdx the first index of the region
	 * @param lastIdx the last locus index of the region
	 */
	int Init(int yOffset, uint firstIdx, uint lastIdx);

	/**
	 * @brief Called for each block to be added to the plot. if we aren't doing labels, we'll render block densities
	 * @param block Haplotype block to be added
	 * @param classification Used to determine how dark the border of the block will be
	 */
	void AddBlock(BlockListNode *block, uint classification);
	//void AddBlock(HaplotypeBlock *block, uint classification);

	/**
	 * @brief Set the density bounds so that the block densities will be correct
	 */
	void SetBlockDensityBounds(uint min, uint max);



protected:
	int curSnpIdx;							///<Keep up with where we are. Assuming that each locus has been ordered
	size_t regionHeight;					///<This is the height of the shaded header
	int upperY;								///<The Y coordinate for the top of the label section
	int lowerY;								///<The bottom for the label section
	float density_range;					///<The range of all block densities 
	int firstSnpIndex;						///<Index of the first snp
	int lastSnpIndex;		
	int minBDensity;						///<The very lowest possible block density expected

};

}

}

#endif
