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

namespace Utility {
bool ConfigurationParser::ParseLine(const char*line, uint linenumber) {
	bool success=false;

	//Skip empty lines or comments
	if (strlen(line)>0 && ((line[0] >= 'A' && line[0] <= 'Z') || (line[0] >= 'a' && line[0] <= 'z'))) {
		string tokenString = line;
		char_separator<char> sep(" \t");
		cfgtokenizertype tok(tokenString, sep);					
		cfgtokenizertype::iterator i = tok.begin();
		cfgtokenizertype::iterator end = tok.end();
		string key = (*i).c_str();
		i++;
		string value = (*i).c_str();
		SetValues(key.c_str(), value.c_str());
		success=true;
	}
	return success;
}


bool ConfigurationParser::GetBoolean(const char *value) {
	bool isOn=false;
	isOn = strcmp(value, "ON")==0 || strcmp(value, "on")==0 || strcmp(value, "On")==0 || strcmp(value, "YES")==0 || strcmp(value, "Yes") ==0 || strcmp(value, "Yes")==0;
	return isOn;
}


}
