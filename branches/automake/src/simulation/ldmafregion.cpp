//
// C++ Implementation: ldmafregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldmafregion.h"

namespace Simulation {

namespace Visualization {

LdMAFRegion::LdMAFRegion(ImageParameters *param) : LdPngComponent(param), curSnpIdx(0), lastXMAF(0), lastYMAF(0) {
	regionHeight = (size_t) (param->headerHeight *0.40);
	regionOffset = param->headerHeight - regionHeight;
}


LdMAFRegion::~LdMAFRegion()
{
}

int LdMAFRegion::Init(int yOffset, uint first, uint last) {
	startY = yOffset;
	stopY  = regionHeight;
	return yOffset - stopY;
}

void LdMAFRegion::AddLocus(Locus& locus) {
	if (imgParams->showMAF) {
		xy pos = imgParams->GetHeaderPosition( curSnpIdx++ );

		size_t curX = pos.x;
		size_t curY = (size_t)(regionHeight - 2.0*locus.GetMinAlleleFreq()*regionHeight);

		if (curSnpIdx -1 % imgParams->snpWidth >  0)
			imgParams->png->line(lastXMAF, pos.y + startY - lastYMAF, curX, pos.y + startY - curY, 0.0,0.0,0.25);
	
		lastXMAF = curX;
		lastYMAF = curY;
	}

}

void LdMAFRegion::Open(const char *name) {
	size_t fontSize = (size_t)(imgParams->fontSize * 0.7);
	size_t x1 = imgParams->margin.x  + imgParams->halfWidth;
	size_t x2 = imgParams->dimensions.x - x1;
	size_t halfFontSize = (uint)(fontSize * 0.5);

	if (imgParams->showMAF) {
		for (size_t i = 0; i<imgParams->wrapCount; i++) {
			xy pos = imgParams->GetHeaderPosition( i*imgParams->snpWidth + 1);
			size_t y1 = pos.y + startY - halfFontSize;
			size_t y2 = y1 - stopY;
			size_t segment = (size_t)(stopY * 0.2);
		

			imgParams->png->line(x1, y1 + halfFontSize, x2, y1 + halfFontSize, 0, 0, 0);
			imgParams->png->line(x1, y2 + halfFontSize, x2, y2 + halfFontSize, 0, 0, 0);

			imgParams->PlotText("MAF", (int)(2+imgParams->fontSize * 1.5), (int) (y2 + (regionHeight*0.5)), 1.570795);

/*			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, 
					(int)(2 + imgParams->fontSize * 1.5), (int)(y2 + (regionHeight*.5)) , 1.570795, "MAF", 0, 0, 0);
	*/
			//int xOffset = (int)(imgParams->margin.x-imgParams->png->get_text_width((char *)imgParams->font.c_str(), fontSize, "0.5")-imgParams->halfWidth * 0.15);
			int xOffset = imgParams->margin.x - 3 - imgParams->png->get_text_width((char *)imgParams->font.c_str(), fontSize, (const char *)"0.5");
//			cout<<"MAF Region: "<<name<<"\t"<<imgParams->margin.x<<"\t"<<
/*			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y1, 0.0, "0.5", 0, 0, 0);
			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y1-segment, 0, "0.4", 0, 0, 0);
			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y1-segment*2, 0, "0.3", 0, 0, 0);
			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y1-segment*3, 0, "0.2", 0, 0, 0);
			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y1-segment*4, 0, "0.1", 0, 0, 0);
			imgParams->png->plot_text((char *)imgParams->font.c_str(), imgParams->fontSize, xOffset, y2, 0, "0.0", 0, 0, 0);
*/	
			imgParams->PlotText("0.5", xOffset, y1);
			imgParams->PlotText("0.4", xOffset, y1-segment);
			imgParams->PlotText("0.3", xOffset, y1-segment*2);
			imgParams->PlotText("0.2", xOffset, y1-segment*3);
			imgParams->PlotText("0.1", xOffset, y1-segment*4);
			imgParams->PlotText("0.0", xOffset, y2);
		
			imgParams->png->filledsquare_blend(x1,y1-segment  + halfFontSize,x2, y1-segment*2+halfFontSize, 0.45, 0.85, 0.85, 0.95);
			imgParams->png->filledsquare_blend(x1,y1-segment*3+halfFontSize,x2, y1-segment*4+halfFontSize, 0.45, 0.85, 0.85, 0.95);
		}

			
	}

}
void LdMAFRegion::Close() {
}


}

}
