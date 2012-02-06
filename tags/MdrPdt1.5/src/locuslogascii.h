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
#ifndef ESELOCUSLOG_H
#define ESELOCUSLOG_H

#include "utility/utility.h"
#include "genetics/locuslog.h"

namespace MDR {
using namespace std;
using namespace Genetics;
using namespace Genetics::Reporting;


/**
Is responsible for opening the locus log and determining what and when to write to it.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LocusLogAscii : public LocusLog, public AsciiLog{
public:
	typedef enum {
		WriteAll,				//All loci reported
		WriteIncludes,			//Only loci that are included in search
		WriteExcludes			//Only loci that have been excluded
	} LocusWriteStyle;

    LocusLogAscii(LocusWriteStyle style, const char *reportname, const char *ext);
    ~LocusLogAscii();

	/**
	 * @brief Log interface
	 */
	
	/**
	 * @brief Writes the snp out to the log
	 */
	void WriteSnp(SnpAligned *snp);

	/**
	 * @brief Returns a basic summary of the configuration details
	 */
	string ReportConfiguration();

	/**
	 * @brief Writes the logs header to the stream
	 */
	void WriteHeader();

	void Close();
	
protected:
	LocusWriteStyle style;				//Used to determine which loci are reported
};


inline
string LocusLogAscii::ReportConfiguration() {
	char cfg[4096];
	string style="WriteAll";

	if (this->style==WriteIncludes)
		style="Included Loci Only";
	else if (this->style==WriteExcludes)
		style="Excluded Loci Only";

	sprintf(cfg, "* Locus Style:             %s\n* Locus Filename:          %s\n", style.c_str(), filename.c_str());
	return cfg;
}

inline
LocusLogAscii::~LocusLogAscii() {
}

inline
LocusLogAscii::LocusLogAscii(LocusWriteStyle style, const char *reportname, const char *ext) 
		: AsciiLog(reportname, ext, 1), style(style) {}

inline
void LocusLogAscii::WriteHeader() {
	*stream<<"Line-Number Locus-Number Code-0/AA Code-1/Aa Code-2/aa Inclusion-State\n";
}

inline
void LocusLogAscii::Close() {
	CloseFile();		
}
}
#endif
