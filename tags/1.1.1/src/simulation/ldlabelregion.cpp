//
// C++ Implementation: ldlabelregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldlabelregion.h"

namespace Simulation {

namespace Visualization {

LdLabelRegion::LdLabelRegion(ImageParameters *param) : LdPngComponent(param), curSnpIdx(0) {
	regionHeight = param->headerHeight;
}


LdLabelRegion::~LdLabelRegion() {
}

int LdLabelRegion::Init(int yOffset, uint firstIdx, uint lastIdx) {
	upperY = yOffset;
	lowerY = regionHeight - yOffset;
	firstSnpIndex = firstIdx;
	lastSnpIndex = lastIdx;

	return lowerY;
}

void LdLabelRegion::AddLocus(Locus& locus) {
	if (imgParams->showLabels) {
		xy pos = imgParams->GetHeaderPosition( curSnpIdx++ );

		//Shade alternate loci differently to set them apart
		if (curSnpIdx%2 > 0)
			//imgParams->png->filledsquare_blend(pos.x-imgParams->halfWidth,pos.y+regionHeight,
			imgParams->png->filledsquare_blend(pos.x-imgParams->halfWidth,pos.y + upperY,
					pos.x+imgParams->halfWidth, pos.y - imgParams->halfWidth, 
					0.85, 0.85, 0.85, 0.95);
		else
			imgParams->png->filledsquare_blend(pos.x-imgParams->halfWidth,pos.y + upperY,
					pos.x+imgParams->halfWidth, pos.y - imgParams->halfWidth, 
					0.85, 0.95, 0.95, 0.95);
		//Write the label to the plot
		imgParams->PlotText(locus.GetLabel().c_str(), (int)(pos.x+(imgParams->fontSize*0.5)), pos.y, 1.570795);
/*		imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, 
					(uint)(pos.x+(imgParams->fontSize * 0.5)), pos.y, 
					1.570795, (char *)locus.GetLabel().c_str(), 
					0.0, 0.0, 0.0);	
*/
	}
}





void LdLabelRegion::SetBlockDensityBounds(uint min, uint max) {
	minBDensity = min;
	density_range = max - min;
	
}


void LdLabelRegion::AddBlock(BlockListNode *block, uint classification) {
//AddBlock(HaplotypeBlock *block, uint classification) {
	//If we aren't doing labels, we will replace them with block densities
	if (!imgParams->showLabels) {
		int nextIdx = block->GetFirstIndex() - firstSnpIndex;
		int lastIdx = block->GetLastIndex() - firstSnpIndex;

/*		cout<<"Attempting to add block to LD Plot: \n";
		cout<<"\tFirst Snp Index    "<<firstSnpIndex         <<"\tLast Snp Index:    "<<lastSnpIndex<<"\n";
		cout<<"\tblock::First Index "<<block->GetFirstIndex()<<"\tblock::Last Index: "<<firstSnpIndex<<"\n";
		cout<<"\tlastIdx:           "<<lastIdx               <<"\tnexIdx:            "<<nextIdx<<"\n";
*/
	
		//Let's get the density relative to the rest of the blocks on the chromosome
		float relDensity = 0.0;
		if (density_range > 0.0) 
			relDensity = (float)(block->GetBlockDensity() - minBDensity) / density_range;

		size_t headerOffsetTop  = (size_t)(regionHeight * 0.5);
		int topY 				= headerOffsetTop + (size_t)(headerOffsetTop * relDensity);

		//And plot the density as a series of blocks over the width of the block (this block might wrap around)
		while (nextIdx >= 0 && nextIdx <= lastIdx ){ //&& nextIdx < lastSnpIndex) {
			xy nextXY = imgParams->GetHeaderPosition( nextIdx++ );
			imgParams->png->filledsquare(nextXY.x - imgParams->halfWidth, nextXY.y + topY, nextXY.x+imgParams->halfWidth, nextXY.y, 0, 0, 0);
		}

	}
}

void LdLabelRegion::Open(const char *name) {
	if (!imgParams->showLabels) {
		size_t headerOffsetTop   = (size_t)(regionHeight * 0.5);

		for (uint i=0; i<imgParams->wrapCount; i++) {
			xy pos = imgParams->GetHeaderPosition( imgParams->snpWidth * i + 1);
			imgParams->png->line(imgParams->margin.x, pos.y + headerOffsetTop, imgParams->dimensions.x - imgParams->margin.x, pos.y + headerOffsetTop, 0.0,0.0,0.25);
		}
	}
}
void LdLabelRegion::Close() {
}



}

}
