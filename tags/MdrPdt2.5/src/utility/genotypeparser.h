//
// C++ Interface: genotypeparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GAWGENOTYPEPARSER_H
#define GAWGENOTYPEPARSER_H

#include <iostream>

namespace Utility {

using namespace std;

/**
Base class for the different parsers for various genotypes

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenotypeParser{
public:
    GenotypeParser();

    virtual ~GenotypeParser();


	/**
	 * Loads the file, filename
	 * @param filename File containing phenotypes
	 * @return true if the file was parsed without fatal errors
	 */
	virtual bool Load(const char* filename) = 0;

	/**
 	 * Dumps a Text Report of the mappings and other details to the stream
	 */
	virtual void GenerateReport(ostream& os, bool complete) = 0;

	/**
	 * Dumps the mapping details to the stream
	 */
	virtual void GenerateMap(ostream &os) = 0;

};

}

#endif
