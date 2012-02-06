//
// C++ Implementation: configurationparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "configurationparser.h"
#include "string.h"

namespace Utility {
bool ConfigurationParser::ParseLine(const char*line, uint linenumber) {
	bool success=false;

	stringstream str(line);
	
	//Skip empty lines or comments
	if (strlen(line)>0 && ((line[0] >= 'A' && line[0] <= 'Z') || (line[0] >= 'a' && line[0] <= 'z'))) {
		string key;
		str>>key;
		string values;
		str>>values;
		SetValues(key.c_str(), values.c_str(), line);
		success=true;
	}
	return success;
}


bool ConfigurationParser::GetBoolean(const char *value) {
	bool isOn=false;
	isOn = strcmp(value, "ON")==0 || strcmp(value, "on")==0 || strcmp(value, "On")==0 || strcmp(value, "YES")==0 || strcmp(value, "Yes") ==0 || strcmp(value, "Yes")==0;
	return isOn;
}

uint ConfigurationParser::CountColumns(const char *line) {
	stringstream ss(line);
	string vals;

	uint count = 0;

	while (!ss.eof()) {
		vals = "xvx";
		ss>>vals;	
		
		if (vals != "xvx")
			count++;
	}
	
	return count;
}



}
