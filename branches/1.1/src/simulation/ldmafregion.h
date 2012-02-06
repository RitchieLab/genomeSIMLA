//
// C++ Interface: ldmafregion
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDMAFREGION_H
#define SIMULATION_VISUALIZATIONLDMAFREGION_H

#include "ldpngcomponent.h"
#include "locus.h"

namespace Simulation {

namespace Visualization {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdMAFRegion : public LdPngComponent {
public:
    LdMAFRegion(ImageParameters *param);

    virtual ~LdMAFRegion();
	//Causes the locus' location be added tot he bar
	void AddLocus(Locus &l);

	void Open(const char *name);
	void Close();

	int Init(int yOffset, uint first, uint last);

protected:
	int curSnpIdx;		
	size_t regionHeight;
	size_t regionOffset;		
	size_t lastXMAF;
	size_t lastYMAF;
	int startY; 
	int stopY;
};

}

}

#endif
