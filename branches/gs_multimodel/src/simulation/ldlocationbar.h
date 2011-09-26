//
// C++ Interface: ldlocationbar
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDLOCATIONBAR_H
#define SIMULATION_VISUALIZATIONLDLOCATIONBAR_H

#include "ldpngcomponent.h"
#include "locus.h"

namespace Simulation {

namespace Visualization {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdLocationBar : public LdPngComponent {
public:
	/**
	 * @brief Constructs the location bar
	 */
	LdLocationBar(ImageParameters *param, size_t locStart, size_t locStop, bool traceToLabel);
	
	/**
	 * @brief Destructor
	 */
	~LdLocationBar();


	//Causes the locus' location be added tot he bar
	void AddLocus(Locus &l);
	void Open(const char *name);
	void Close();
	int Init(int yOffset, uint first, uint last);
	void SetTracerToLabel(bool doTrace) { traceToLabel=doTrace; }


protected:
	int curSnpIdx;					///<Keep up with where we are. Assuming that each locus has been ordered
	int markerLength;				///<How long the location marker line is
	bool traceToLabel;				///<Determines if we want tracers drawn
	uint segmentStart;				///<Starting location of the segment
	uint segmentStop;
	float segmentWidth;				///<How long this segment is (Megabases)	


	int y1;						///<This is the top of the bar
	int y2;						///<The bottom of the bar
	int y3;						///<This is the actual position of the label box

//	uint headerOffset;				///<The distance from the top of the page

};

}

}

#endif
