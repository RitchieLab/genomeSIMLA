//
// C++ Interface: matchedoddsratio
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEMATCHEDODDSRATIO_H
#define ESEMATCHEDODDSRATIO_H
#include "utility/utility.h"
#include <ostream>

namespace MDR {

using namespace std;

/**
@brief Perform the matched odds ratio

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class MatchedOddsRatio{
public:
    MatchedOddsRatio();

    ~MatchedOddsRatio();

	/**
	 * @brief Calculates the appropriate pairs based on the observed values from a single family and adds them to the local counts
	 * @param truePos  - HR Affected
	 * @param falsePos - LR Affected
	 * @param falseNeg - HR Unaffected
	 * @param trueNeg  - LR Unaffected
	 */
	void AddFamily(uint truePos, uint falsePos, uint falseNeg, uint trueNeg);
	
	/**
	 * @brief Write the ratio and it's table to the report
	 */
	void GenerateReport(ostream *os);

	/**
	 * @brief Sets the details back to 0
	 */
	void Reset();

	/**
	 * @brief toggle verbose mode on and off
	 */
	static bool Verbose;
	
	float GetRatio();
protected:
	uint exposures[4];				//Just easy way to keep up with the data
	uint familyGroups;
	stringstream setupDetails;
};

inline
MatchedOddsRatio::MatchedOddsRatio() {
	Reset();
}

inline
MatchedOddsRatio::~MatchedOddsRatio() { }

inline
void MatchedOddsRatio::Reset() {
	familyGroups = 0;
	for (uint i=0; i<4; i++) 
		exposures[i]=0;
}

}

#endif
