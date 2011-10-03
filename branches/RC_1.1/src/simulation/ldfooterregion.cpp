//
// C++ Implementation: ldfooterregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldfooterregion.h"

namespace Simulation {

namespace Visualization {



LdFooterRegion::LdFooterRegion(ImageParameters *param) : LdPngComponent(param) {
	fontSize = param->labelFontSize;
	pngWriter = param->png;
	font = param->font;
	label = "Chromosome Z: No Label Given";
}


LdFooterRegion::~LdFooterRegion() {
}

int LdFooterRegion::Init(int yOffset, uint firstIdx, uint lastIdx) {
	return imgParams->dimensions.y;
}

void LdFooterRegion::Close() {
	int fontSize = imgParams->labelFontSize;
	int textWidth = pngWriter->get_text_width((char*)font.c_str(), fontSize, (char*)label.c_str());
	int left = (imgParams->dimensions.x / 2) - (textWidth / 2);

	while (left < 0) {
		fontSize = (int)(fontSize * 0.9);
		int textWidth = pngWriter->get_text_width((char*)font.c_str(), fontSize, (char*)label.c_str());
		left = (imgParams->dimensions.x / 2) - (textWidth / 2);
	}

	imgParams->PlotText(label.c_str(), left, (int)(fontSize*0.50));
/*	imgParams->png->plot_text((char*)font.c_str(), fontSize, 
				left, (int)(fontSize * 0.50), 0.0, (char*)label.c_str(), 0.0, 0.0, 0.0);	
 */

}


}

}
