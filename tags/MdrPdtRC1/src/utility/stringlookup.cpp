//
// C++ Implementation: stringlookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "stringlookup.h"

namespace Utility {

StringLookup::~StringLookup()
{
}

uint StringLookup::GetValue(const char *key) {
	uint val = 0;	
	if (key) {
		val = lookupTable[key];
		if (val == 0)	{
			lookupTable[key] = val = nextVal++;
			
		}
	}
	return val;
}

void StringLookup::Report(ostream& os)	{
	os<<"\n***** "<<headerName<<" *****\n";
	
	LookupType::iterator end = lookupTable.end();
	
	for (LookupType::iterator itr = lookupTable.begin(); itr!=end; itr++) 
		os<<"\t"<<itr->first<<" -> "<<itr->second<<"\n";

	os<<"\n";
}

uint StringLookup::Size() {
	return nextVal - 1;
}


StringLookup::StringLookup() : nextVal(1)
{
}

void StringLookup::SetHeader(const char* name) {
	headerName = name;
}

}
