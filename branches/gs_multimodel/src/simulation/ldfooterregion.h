//
// C++ Interface: ldfooterregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDFOOTERREGION_H
#define SIMULATION_VISUALIZATIONLDFOOTERREGION_H
#include "ldpngcomponent.h"

namespace Simulation {

namespace Visualization {

/**
@brief Writes some general text to the image at the bottom of the picture

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdFooterRegion : public LdPngComponent {
public:
	LdFooterRegion(ImageParameters *param);
	virtual ~LdFooterRegion();

	int Init(int yOffset, uint firstIndex, uint lastIdx);
	void Close();
	void Open(const char *) { }
	
	void SetLabel(const char *label) { this->label=label; }
protected:
	string label;
	pngwriter *pngWriter;
	string font;
	int fontSize;
};

}

}

#endif
