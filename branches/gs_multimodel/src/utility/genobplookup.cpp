//
// C++ Implementation: genobplookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genobplookup.h"

#include <assert.h>

namespace Utility {



/**
 * Returns the mapped value setting it if necessary
 */
int GenoBPLookup::GetValue(const char* key)	{
	int value = 0;
	char *str = (char *)key;

	//assert(strlen(key) == 2);
	
	if (str[0] == str[1])
		value = GetHomoZygote(key);
	else
		value = GetHeteroZygote(key);

	return value;
}



int GenoBPLookup::GetHomoZygote(const char* value) {
	int returnValue = NODATA;

	if (strncmp(homozygote1.c_str(), value, 2) == 0)
		returnValue = HOMOZYGOTE1;
	else if (strncmp(homozygote2.c_str(), value, 2) == 0)
		returnValue = HOMOZYGOTE2;
	else if (homozygote1 == "" && (value[0]=='T' || value[0]=='G'))	{
		returnValue = HOMOZYGOTE1;
		homozygote1 = value;
	} else if (homozygote2 == "") {
		returnValue = HOMOZYGOTE2;
		homozygote2 = value;
	}
	else if (homozygote1=="") {
		returnValue = HOMOZYGOTE1;
		homozygote1 = value;
	} 		
	else	{
		nodata.push_back(value);
		printf("Skipping genotype: %s\n", value);
	}
	return returnValue;
}

}
