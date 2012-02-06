//
// C++ Interface: basiclog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYBASICLOG_H
#define UTILITYBASICLOG_H
#include <cstdarg>
#include <fstream>
#include <iostream>
#include "types.h"

namespace Utility {

using namespace std;

/**
@brief Base class for simple logs

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BasicLog{
public:
    BasicLog();
    virtual ~BasicLog();
	
	/**
	 * @brief Write bytes to the stream
	 */
	virtual void Write(const char *data) = 0;

	/**
	 * Writes formatted data to the file
	 */
	virtual void Write(uint level, const char *format, ...) = 0;
	
	/**
	 * @brief Return a pointer to the stream. This can be COUT depending on the log type and how it's set up
	 */
	virtual ostream *GetStream() = 0;

	/**
	 * @brief Returns the "filename" associated with the log
	 */
	virtual string GetFilename() = 0;

protected:
};

class AsciiLog : public BasicLog {
public:
	/**
	 * @brief Set up a log with a report name and extension
	 */
	AsciiLog(const char *report, const char *ext, uint level);

	/**
	 * @brief set up the report when you know the filename
	 */
	AsciiLog(const char *filename, uint level);
	~AsciiLog();

	/**
	 * @brief Forces a write to the log
	 */
	void Write(const char *data);
	/**
	 * @brief Writes formatted data to the file
	 */
	void Write(uint level, const char *format, ...);

	/**
	 * @brief Returns the stream associated with the log
	 */
	ostream *GetStream();

	/**
	 * @brief Opens the stream up and sets up necessary paramenters. Filename is assumed to be a valid filename
	 */
	ostream *Open(const char *filename);

	virtual void CloseFile();

	/**
	 * @brief Returns the name of the file associated with the log
	 */
	string GetFilename();
protected:
	uint level;					///<Can be used to control the level of output
	string filename;			///<Name of file being written to
	ofstream *file;				///<File being written to
	ostream *stream;			///<Stream being written to
};


inline
string AsciiLog::GetFilename() {	
	if (strlen(filename.c_str())==0)
		return "Standard Out\n";
	return filename;
}

inline
ostream *AsciiLog::GetStream() {
	return stream;
}

inline
void AsciiLog::Write(const char *data) {
	if (stream)
		*stream<<data;
} 

inline
void AsciiLog::Write(uint level, const char *format, ...) {
	if (level > this->level && stream)
		return;

	va_list args; 						///<iterate over the arguments*/
	int integer;						///<used for integer args
	double dbl;							///<Used for doubles
	char *p=(char *)format;				///<
	va_start(args,format); /* point to first element after fmt*/
	while(*p)	{
		if (strncpy(p, "%%", 2) == 0) {
			++p;
			if (*p=='i') {
				integer=va_arg(args,int);
				*stream<<integer;
			}
			else if (*p=='d') 	{
				dbl=va_arg(args,double);
				*stream<<dbl;
			}	
			else	{
				stream->write(p, 1);
			}
		}
		else {
			stream->write(p, 1);
		}	
				
		++p; /* get the next char of the format string */
	}/*while*/
	va_end(args); /*cleanup*/ 	
}


inline
ostream *AsciiLog::Open(const char *filename) {
	file = new ofstream(filename,ios::out|ios::trunc);
	stream = file;
	if (!file->good()) {
		delete file;
		file = NULL;	
		stream = &cout;
		*stream<<"Unable to open file: "<<filename<<". Using standard out.\n";
	}
	else 
		this->filename = filename;
	return stream;

}

inline
AsciiLog::AsciiLog(const char *filename, const char *ext, uint level) : level(level), stream(NULL) {
	if (strlen(ext) == 0 || strcmp(ext, "NOLOG") == 0) {
		stream = NULL;
		file = NULL;
		this->filename = "No Report Defined";
	}
	else if (strcmp(ext, "STDOUT") == 0 || strlen(filename) == 0) {
		stream = &cout;
		file = NULL;
		this->filename = "Standard Out";
	}
	else {
		this->filename = string(filename) + string(".") + string(ext);
		Open(this->filename.c_str());
	}	
}

inline
AsciiLog::AsciiLog(const char *filename, uint level) : level(level), filename(filename), stream(NULL) {
	if (strcmp(filename, "STDOUT")==0) {
		stream = &cout;
		file = NULL;
		this->filename = "Standard Out";		
 	}
	else if (strcmp(filename, "NOLOG") == 0) {
		stream = NULL;
		file = NULL;
		this->filename = "No Report Defined";
	}
	else 
		Open(filename);
} 

inline
void AsciiLog::CloseFile() {
	if (file) {
		file->close();
		delete file;
		file = NULL;
	}
}
		
inline
AsciiLog::~AsciiLog() {
	if (file) {
		file->close();
		delete file;
		file = NULL;
	}
}		
		

		
		
}

#endif
