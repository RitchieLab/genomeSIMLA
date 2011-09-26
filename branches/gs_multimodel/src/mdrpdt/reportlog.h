//
// C++ Interface: reportlog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEREPORTLOG_H
#define ESEREPORTLOG_H

#include <string>
#include <ostream>
#include "genetics/snpcontainer.h"

namespace MDR {

using namespace std;
using namespace Genetics;
/**
@brief Base class for the different report files. 
This class has one functional role other than serving as an interface guide. ReportLog::Create(ReportStyle) will return the 
appropriate object type.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ReportLog{
public:

    virtual ~ReportLog();

	typedef enum {
		NoReport=0,									//0 - No Reporting
		PerfectOnly,								//1 - Perfect Models Only
		AboveThresholdNoMetric,						//2 - Above MD Threshold but without the metric
		AboveThresholdWithMetric,					//3 - Above MD Threshold including metric
		TheNumberFour,								//4 - Above MD Threshold including metric including some additional statistical information
		TheNumberFive,								//5 - I have no clue! This isn't well documented
		Binary=6 									//6 - In binary format
	} ReportStyle;

	/** 
	 * @brief Factory functionality
	 */
	static ReportLog *Create(const char *type, ReportStyle style, const char *filename);
	/**
	 * Report Interface
	 */

	/**
	 * @brief Returns a pointer to the stream
	 * I don't really like this approach, but there are a few things that need it this way right now
	 */
	virtual ostream *GetStream()=0;
	
	/**
	 * @brief Initialize the log. This will open the file and set up any other general details
	 * @return Error Message or "" if no error is encountered
	 */
	virtual string Init() = 0;

	/**
	 * @brief Generates string to display configuration for the user
	 */
	virtual string ReportConfiguration() = 0;

	void ReportOnRepository(const char *title, uint locusCount, SnpReporter *snps);
	
	virtual void EchoTitle(const char *title)=0;

	//Write plain text to the report
	virtual void Echo(const char *entry)=0;

	//Write the text to the report accompanying a value
	virtual void Echo(uint val, const char *entry)=0;
	
	virtual void WriteHeader(uint locusCount)=0;
	virtual void ReportNoModels(const char *)=0;

protected:
	string reportType;
	/**
	 * Inaccessible interface
	 */
    ReportLog(const char *type);
};


inline
void ReportLog::ReportOnRepository(const char *title, uint locusCount, SnpReporter *snps) {
	char entryText[516];
	for (uint loci=0; loci<locusCount; loci++) {
		if (locusCount>0) 
			sprintf(entryText, "%d SNP-models %d/%d models reported/found meeting or exceeding threshold", loci+1, snps->GetEntryCount(loci, 0), snps->GetTotalEntryCount(loci, 0));
		else  
			sprintf(entryText, "%d/%d models reported/found meeting or exceeding threshold", snps->GetEntryCount(loci, 0), snps->GetTotalEntryCount(loci, 0));
		
			Echo(entryText);
		uint count=snps->GetEntryCount(loci, 0);
		if (count > 0) {
			EchoTitle(title);
			WriteHeader(loci);
			for (uint i=0; i<count; i++) {
				Echo(i, snps->GetReportEntry(0, loci, i).c_str());
			}
		}
		else {
			ReportNoModels(title);
		}
	}
}


}

#endif
