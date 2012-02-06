//
// C++ Interface: locuslog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_REPORTINGLOCUSLOG_H
#define GENETICS_REPORTINGLOCUSLOG_H

#include "snpaligned.h"

namespace Genetics {

namespace Reporting {



/**
@brief Base class for locus reports

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LocusLog {
public:
	LocusLog() {}
	virtual ~LocusLog() {}
	/**
	 * @brief Writes the snp out to the log
	 */
	virtual void WriteSnp(SnpAligned *snp) = 0;

	virtual string ReportConfiguration() =0;

	virtual void WriteHeader() = 0;

	virtual void Close() =0;
};

}

}

#endif
