//
// C++ Interface: distmappingfn
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONDISTMAPPINGFN_H
#define SIMULATIONDISTMAPPINGFN_H

#include "utility/types.h"
#include <iostream>
#include <math.h>

namespace Simulation {

using namespace std;

/**
Base class for the different distance mapping functions

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class DistMappingFn{
public:
	uint position;

    DistMappingFn() : position(0) { }
    virtual ~DistMappingFn() {}

	void Reset(int pos=0) {
		position = pos;
	}

	virtual string GetType() = 0;

	virtual uint GetMapPosition(float recombFraction)=0;
	virtual void Report(ostream& os) = 0;
};


/**
 * @brief Convert probability to distance based on Haldane function
 * @param recProb Probability of recombination
 */
class HaldaneMapping : public DistMappingFn {
public:
	HaldaneMapping() {}
	~HaldaneMapping() {}

	uint GetMapPosition(float recProb) {
		float Cab2 = (float)(1.0 - (2.0 * recProb));
		double Xab = -0.5 * logf(Cab2);
		position += (uint)(Xab * 100000000.0);
		return position;
	}

	void Report(ostream& os) {
		os<<"Haldane Mapping\n";
	}

	string GetType() { return "HALDANE"; }
};

/**
 * @brief Convert probability to distance based on kosamib equation
 * @param recProb Probability of recombination
 */
class KosambiMapping : public DistMappingFn {
public:
	KosambiMapping() {}
	~KosambiMapping() {}
	uint GetMapPosition(float recProb) {
		float Cab2 = (float)2.0 * recProb;
		double x = ((1.0 + Cab2)/(1.0-Cab2));
		double Xab = (float)0.25 * logf(x);
		position += (uint)(Xab * 100000000.0);
		return position;
	}

	void Report(ostream& os) {
		os<<"Kosambi Mapping\n";
	}

	string GetType() { return "KOSAMBI"; }

};

}

#endif
