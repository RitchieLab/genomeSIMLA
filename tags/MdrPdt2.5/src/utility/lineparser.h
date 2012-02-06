//
// C++ Interface: lineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESELINEPARSER_H
#define ESELINEPARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "types.h"

namespace Utility {

using namespace std;


/**
 * @brief Baseclass for the ascii parser routines
 */
class AsciiParser {
public:
	AsciiParser(){}
	virtual ~AsciiParser() {}
	virtual bool ParseLine(const char *line, uint val)=0;
};


/**
  @brief Breaks the line into a vector and passes it to the functor
 */
class ArrayedTextParser : public AsciiParser {
public:
	ArrayedTextParser() {}
	virtual ~ArrayedTextParser() {} 
	bool ParseLine(const char *line, uint val);
	virtual bool Evaluate(StringArray& words, uint val)=0;
};



/**
@brief Parses ascii files and executes a functor for each line

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LineParser{
public:
    LineParser() {}

    ~LineParser(){}

	/**
	 * @brief Opens a file and sends each line to the function pointer at fn
	 */
	uint Parse(const char* filename, AsciiParser *parser);

};

/**
  @brief Converts a file to an array of strings
 */
class FileToArray : public AsciiParser {
public:
	FileToArray();
	~FileToArray();
	
	bool ParseLine(const char *line, uint val=0);
	StringArray strings;
};


inline
FileToArray::FileToArray() {}

inline
FileToArray::~FileToArray() {}

inline
bool FileToArray::ParseLine(const char *line, uint val) {
	strings.push_back(string(line));
	return true;
}

/**
 * @brief A simple parser object that returns the number of lines encountered. 
 * Eventually, this could be expanded to get other types of information like the number of columns in certain types of lines
 */
class BasicLineCounter: public AsciiParser {
public:
	BasicLineCounter() {}
	~BasicLineCounter() {}

	bool ParseLine(const char *line, uint val=0);
};

inline
bool BasicLineCounter::ParseLine(const char *line, uint val) {
	if (strlen(line) > 0 && line[0] != '/' && line[0] != '#')
		return true;
	else
		return false;
}

/**
@brief Sets values in the bitvector to value if they appear in the file
*/
class ExclusionList : public AsciiParser {
public:
	ExclusionList(BitSetType *active, bool value=false) : active(active), value(value) {}
	~ExclusionList() {}
	
	bool ParseLine(const char *line, uint val=0);
protected:
	BitSetType *active;
	bool value;

};


inline
bool ExclusionList::ParseLine(const char* line, uint val /*=0*/) {
	bool success = strlen(line) > 0;

	if (success)	{
		uint idx=atoi(line);	
		(*active)[idx-1]=value;
	}
	return success;
}




inline
bool ArrayedTextParser::ParseLine(const char *rawline, uint val) {
	StringArray line; 
	uint len=strlen(rawline);
	
	if (len > 0) {
		string s=rawline;
	
		boost::char_separator<char> sep(" \t,\r", "", boost::drop_empty_tokens);;

		strtokenizer tok(s, sep);						
		strtokenizer::iterator i = tok.begin();
		strtokenizer::iterator end = tok.end();
		
		do  {
			line.push_back(*i);
		}
		while (++i != end) ;
			
		return Evaluate(line, val);
	}	
	return false;
}








}

#endif
