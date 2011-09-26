//
// C++ Interface: reportlogbinary
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEREPORTLOGBINARY_H
#define ESEREPORTLOGBINARY_H

#include "reportlog.h"
#include <assert.h>
#include <iostream>
#include <fstream>


namespace MDR {

using namespace std;

/**
@brief Binary report file

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ReportLogBinary : public ReportLog {
friend class ReportLog;
public:
 	string Init();
	ostream *GetStream();
	string ReportConfiguration();
	void Echo(uint val, const char *entry);
	void Echo(const char *entry);
	void EchoSnp(uint val, SnpAligned *snp);
	void EchoTitle( const char *title);
	void WriteHeader(uint locusCount);
	void ReportNoModels( const char *title);
    ~ReportLogBinary();
protected:
    ReportLogBinary(const char *type, ReportLog::ReportStyle style, const char *filename);
	ReportLog::ReportStyle style;			//Right now, this is pretty much meaningless since there is one setting for binary
	ostream *log;
	const char *filename;
};

inline
ReportLogBinary::ReportLogBinary(const char *type, ReportLog::ReportStyle style, const char *filename) : ReportLog(type), style(style), log(NULL), filename(filename) {}

inline
ReportLogBinary::~ReportLogBinary() {
	if (style!=ReportLog::NoReport) {
		((ofstream*)log)->close();	
		delete log;	
	}
}

inline
void ReportLogBinary::ReportNoModels( const char *title) {
	cerr<<"BINARY REPORT NOT SETUP\n";
	assert(0);
}

inline
void ReportLogBinary::Echo(uint val, const char *entry) {
	cerr<<"BINARY ECHO NOT SETUP YET\n";
}
inline
void ReportLogBinary::Echo(const char *entry) {
	cerr<<"BINARY ECHO NOT SETUP YET\n";
}
inline
void ReportLogBinary::EchoSnp(uint val, SnpAligned *snp) {
	cerr<<"BINARY ECHO NOT SETUP YET\n";
	assert(0);
}	

inline
void ReportLogBinary::EchoTitle( const char *title) {
	//This doesn't really mean much for a highly compressed report
}

inline
void ReportLogBinary::WriteHeader(uint locusCount) {
	*log<<"TBD";
}

inline
string ReportLogBinary::ReportConfiguration() {
	char cfg[4096];
	string style="Undefined Report Mechanism";
	if (this->style==Binary)
		style="Binary";
		
	sprintf(cfg, "Report Style:            %s\nReport Filename:         %s\n", style.c_str(), filename);
	return cfg;
}

inline
ostream *ReportLogBinary::GetStream() {
	if (log)
		return log;
	else	{
		cerr<<Init();
		assert(log);
		return &cout;
	}
}

inline
string ReportLogBinary::Init() {
	string errMsg = "";
	if (style == ReportLog::NoReport)
		log = &cout;
	else {
		log=new ofstream(filename,ios::out|ios::trunc|ios::binary);
		if (!((ofstream*)log)->good())
			errMsg="There was a problem opening the file: "+string(filename);
	}
	return errMsg;
}

}

#endif
