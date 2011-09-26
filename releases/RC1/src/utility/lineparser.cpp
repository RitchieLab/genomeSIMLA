//
// C++ Implementation: lineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "lineparser.h"

namespace Utility {

uint LineParser::Parse(const char *filename, AsciiParser *parser) {
	char line[MAX_LINE_LENGTH];
	uint linesParsed = 0;							//Keep up with what we've seen
	uint lineCount	= 0;							//Used to know how many were skipped

	ifstream file(filename, ios_base::in);			//open file for reading
	if (file.is_open()) {
		while (!file.eof()) {
			file.getline(line, MAX_LINE_LENGTH);
			if (parser->ParseLine(line, linesParsed))
				lineCount++;
			linesParsed++;
		}
	}
	else
		cout<<"File "<<filename<<" wasn't opened\n";
	return lineCount;
}
}//Utility
