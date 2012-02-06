//
// C++ Interface: gtfileparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_PARSERGTFILEPARSER_H
#define GENETICS_PARSERGTFILEPARSER_H

#include "snpaligned.h"

namespace Genetics {

namespace Parser {

/**
@brief Base class for the various genotype File Parsers
It is expected that the one that instantiates the parser knows exactly how to construct it. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GtFileParser{
public:
    GtFileParser();

    virtual ~GtFileParser();

	/**
	 * @brief Returns the next snp
	 * @return SnpAligned pointer
	 */
	virtual SnpAligned *NextSnp()=0;
	/**
	 * @brief Opens up each file in the current group
 	 * @return True/False depending on wither there were problems opening the file(s) (or there were no more to be opened)
	 */
	virtual bool Open()=0;
	/**
	 * @brief Close the currently opened files
	 */
	virtual void Close()=0;

	virtual void Reset()=0;

	/**
	 * @brief Returns the number of snps found in the file(s)
	 */
	virtual uint GetSnpCount()=0;

	/**
	 * @brief Returns the number of individuals found in the file(s)
	 */
	virtual uint GetIndividualCount()=0;

	//virtual uint GetAffectedCount()=0;
	//virtual uint GetUnaffectedCount()=0;
	virtual void GetStatusMask(CaseControlStatus *)=0;			
	virtual void AppendScrambledStatus(CaseControlStatus *)=0;

	virtual void GenerateReport(ostream &os)=0;
};


inline
GtFileParser::GtFileParser() {}

inline
GtFileParser::~GtFileParser() {	
}


}

}

#endif
