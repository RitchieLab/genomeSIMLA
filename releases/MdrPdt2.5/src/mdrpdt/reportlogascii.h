//
// C++ Interface: reportlogascii
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEREPORTLOGASCII_H
#define ESEREPORTLOGASCII_H

#include "reportlog.h"
#include <assert.h>
#include <iostream>
#include <fstream>

namespace MDR {

using namespace std;

/**
ASCII based report file. 
	typedef enum {
		NoReport=0,									//0 - No Reporting
		PerfectOnly,								//1 - Perfect Models Only
		AboveThresholdNoMetric,						//2 - Above MD Threshold but without the metric
		AboveThresholdWithMetric,					//3 - Above MD Threshold including metric
		TheNumberFour,								//4 - Above MD Threshold including metric including some additional statistical information
		TheNumberFive,								//5 - I have no clue! This isn't well documented
		Binary=6 									//6 - In binary format
	} ReportStyle;

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ReportLogAscii : public ReportLog {
public:
    ReportLogAscii(const char *type, ReportLog::ReportStyle style, const char *filename);

	string Init();
	ostream *GetStream();
	string ReportConfiguration();
	void EchoSnp(uint val, SnpAligned *snp);
	void Echo(const char *entry);
	void Echo(uint val, const char *entry);
	void WriteHeader(uint locusCount);
	void EchoTitle(const char *title);

    ~ReportLogAscii();
	void ReportNoModels(const char *title);
protected:
	ReportLog::ReportStyle style;
	ostream *log;
	const char *filename;
};

inline
void ReportLogAscii::WriteHeader(uint locusCount) {
	char metric[20]="";
	char header[1024];

	if (style==AboveThresholdWithMetric)
		strcpy(metric, "MD_Value ");


	sprintf(header, "Ranking Model %sModelID", metric);
	
	*log<<header<<"\n";
}

inline
void ReportLogAscii::ReportNoModels( const char *title) {
	*log<<title<<": \n"<<"None found.\n";
}

inline
void ReportLogAscii::EchoTitle( const char *title) {
	*log<<title<<":\n";
}

inline
void ReportLogAscii::Echo(const char *entry) {
	*log<<entry<<"\n";
}

inline
void ReportLogAscii::Echo(uint val, const char *entry) {
	*log<<val<<" "<<entry<<"\n";
}	


inline
void ReportLogAscii::EchoSnp(uint val, SnpAligned *snp) {
	*log<<val<<" "<<snp->GetLineNumber()<<" ";
	if (style==AboveThresholdWithMetric)
		*log<<snp->GetEvaluation( 0 );
		//*log<<snp->GetLastMdEval()<<" ";
	*log<<snp->GetTxtDescriptor()<<"\n";
}	

inline
string ReportLogAscii::ReportConfiguration() {
	char cfg[4096];
	string style="Undefined Report Mechanism";
	string filename=this->filename;

	if (filename=="")
		filename="Standard Out";

	if (this->style==NoReport)
		style="No Report";
	if (this->style==PerfectOnly)
		style="Only Perfect Models";
	else if (this->style==AboveThresholdNoMetric)
		style="Above Metric Threshold- No Metric Written";
	else if (this->style==AboveThresholdWithMetric)
		style="Above Metric Threshold- Metric Included";
	else if (this->style==Binary)
		style="Binary";
		
	sprintf(cfg, "%s:          %s\nReport Filename:         %s\n", reportType.c_str(), style.c_str(), filename.c_str());
	return cfg;
}


inline
ReportLogAscii::~ReportLogAscii() {
	if (style!=ReportLog::NoReport && strcmp(filename, "")!=0) {
		((ofstream*)log)->close();
		delete log;
	}
}


inline
ReportLogAscii::ReportLogAscii(const char *type, ReportLog::ReportStyle style, const char *filename) : ReportLog(type), style(style), log(NULL), filename(filename) {}

inline
ostream *ReportLogAscii::GetStream() {
	if (log)
		return log;
	else	{
		cerr<<Init();
		assert(log);
		return &cout;
	};
}

inline
string ReportLogAscii::Init() {
	string errMsg = "";
	if (style == ReportLog::NoReport || strcmp(filename,"")==0)
		log = &cout;
	else {
		log = new ofstream(filename,ios::out|ios::trunc);
		if (!((ifstream*)log)->good())
			errMsg="There was a problem opening the file: "+string(filename);
	}
	return errMsg;
}
		

}

#endif
