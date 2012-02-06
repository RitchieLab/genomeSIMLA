//
// C++ Interface: snpcontainer
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESESNPCONTAINER_H
#define ESESNPCONTAINER_H
#include <string>
#include "utility/utility.h"
#include "snpaligned.h"

namespace Genetics {

using namespace std;
/**
 *  @brief Base class for classes that hold report details.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>

 */
class SnpReporter {
public:
	SnpReporter();
	virtual ~SnpReporter();
	
	/**
	 * @brief Returns the number of items logged as opposed to the number of items to be reported
	 */
	virtual uint GetTotalEntryCount(uint modelSize, uint foldCount)=0;

	/**
	 * @brief Returns the number of items to be reported
	 */
	virtual uint GetEntryCount(uint modelSize, uint foldCount)=0;

	/**
	 * @brief Returns the text entry for a given loci/index
	 */
	virtual string GetReportEntry(uint foldID, uint modelSize, uint snpID)=0;
};
/**
@brief Base class for a object to hold snps.
In addition to receiving snps, classes of this type can be used to return them as well.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpContainer : public SnpReporter {
public:
    SnpContainer();

    virtual ~SnpContainer();

	virtual SnpAligned *GetSnp(uint n)=0;
	

};



inline
SnpReporter::SnpReporter() {}

inline
SnpReporter::~SnpReporter() {}

inline
SnpContainer::SnpContainer() {}

inline
SnpContainer::~SnpContainer() {}

}

#endif
