//
// C++ Implementation: ldlocationbar
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldlocationbar.h"
#include <iomanip>

namespace Simulation {

namespace Visualization {

LdLocationBar::LdLocationBar(ImageParameters *param, size_t locStart, size_t locStop, bool traceToLabel) : 
		LdPngComponent(param), curSnpIdx(0), markerLength(5), traceToLabel(traceToLabel), segmentStart(locStart) {
	segmentStop = locStop;

	//segmentWidth = (float)(locStop-locStart)/(float)(imgParams->imageWidth-2*imgParams->halfWidth-2*imgParams->margin.x);	
}


LdLocationBar::~LdLocationBar() {
}

int LdLocationBar::Init(int offsetY, uint first, uint last) {
	segmentWidth = (float)(imgParams->dimensions.x-4*imgParams->halfWidth-2*imgParams->margin.x)/(float)(segmentStop-segmentStart);


	y1 = offsetY;
	y3 = y2 = offsetY - markerLength;
	if (traceToLabel) 
		y3 = y2 - 50;
	return y3;
}


void LdLocationBar::AddLocus(Locus& locus) {
	if (imgParams->showPositions) {
		xy pos = imgParams->GetHeaderPosition( curSnpIdx++ );
		uint origX = (int)(((float)(locus.GetLocation() - segmentStart) * segmentWidth) 
					+ imgParams->margin.x + imgParams->halfWidth * 2);
		size_t y1 = pos.y + this->y1;
		size_t y2 = pos.y + this->y2;
		size_t y3 = pos.y + this->y3;

		//Draw location bar
		imgParams->png->line(origX, y1, origX, y2, 0, 0, 0);

		//imgParams->halfWidth
		if (traceToLabel) 
			imgParams->png->line(origX, y2, pos.x, y3, 0,0,0);
	}
}

void LdLocationBar::Open(const char *name) {
	if (imgParams->showPositions) {
		
		size_t x1 = imgParams->margin.x + imgParams->halfWidth * 2;
		size_t x2 = imgParams->dimensions.x - x1;
	
		for (size_t i = 0; i<imgParams->wrapCount; i++) {
			xy pos = imgParams->GetHeaderPosition( i*imgParams->snpWidth + 1);
			size_t y1 = pos.y + this->y1;
			size_t y2 = pos.y + this->y2;


			imgParams->png->line(x1, y1, x2, y1, 0, 0, 0);
			imgParams->png->line(x1, y2, x2, y2, 0, 0, 0);
		}
	}
}
void LdLocationBar::Close() {
}

}

}
