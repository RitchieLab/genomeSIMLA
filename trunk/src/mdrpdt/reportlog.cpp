//
// C++ Implementation: reportlog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "reportlog.h"
#include "reportlogascii.h"
#include "reportlogbinary.h"

namespace MDR {

ReportLog::ReportLog(const char *type) : reportType(type)
{
}


ReportLog::~ReportLog()
{
}


ReportLog *ReportLog::Create(const char *type, ReportStyle style, const char *filename) {
	if (style < Binary)
		return new ReportLogAscii(type, style, filename);
	else
		return new ReportLogBinary(type, style, filename);
}

}
