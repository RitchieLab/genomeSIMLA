//
// C++ Interface: exception
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYEXCEPTION_H
#define UTILITYEXCEPTION_H

#include <string>
#include <iostream>

namespace Utility {

using namespace std;
/**
Base class for all exceptions

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Exception{
public:
    Exception() {}

    virtual ~Exception() {}

	virtual string GetErrorMessage() const = 0;
	friend ostream & operator << (ostream& os, const Exception& e);
private:

};

inline
ostream &operator<< (ostream& os, const Exception& e) {
	os<<e.GetErrorMessage();
	return os;
}

class FileIOError{
public:
	typedef enum ErrorType {
		FileNotFound=0,
		UnableToRead, 
		UnableToWrite
	} ErrorType;


	FileIOError(const char *filename, FileIOError::ErrorType &type);
	virtual ~FileIOError() {}

	virtual string GetErrorMessage() const;  
protected:
	ErrorType type;
	string filename;
};

inline
FileIOError::FileIOError(const char *filename, FileIOError::ErrorType &type) : type(type), filename(filename) {}

inline
string FileIOError::GetErrorMessage() const {
	stringstream ss;
	switch(type) {
		case FileNotFound:
			ss<<"Unable to find the file "<<filename<<".";
			break;
		case UnableToRead:
			ss<<"Unable to read file, "<<filename<<".";
			break;
		case UnableToWrite:
			ss<<"Unable to write file, "<<filename<<".";
		default:
			ss<<"Unexpected error code: "<<type<<" ("<<filename<<")\n";
	}
	return ss.str();
}

}
#endif
