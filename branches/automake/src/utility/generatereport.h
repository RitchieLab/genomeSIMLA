//
// C++ Interface: generatereport
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYGENERATEREPORT_H
#define UTILITYGENERATEREPORT_H

#include <iostream>

namespace Utility {

using namespace std;

/**
This functor can be used on an STL container to generate a report for each item in the container.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenerateReport{
public:
	/**
	 * Indicates which stream and how you want them separated
	 */
    GenerateReport(ostream &os, const char *del);

    ~GenerateReport();

	template<class T>
	void operator() (T &item) {
		if (printedAlready)
			os<<delim;
		os<<item;
		printedAlready=true;
		//item>>os>>delim;
	}

protected:
	ostream &os;
	string delim;
	bool printedAlready;
};

inline
GenerateReport::GenerateReport(ostream& os, const char *del) : os(os), delim(del), printedAlready(false) {}
}

#endif
