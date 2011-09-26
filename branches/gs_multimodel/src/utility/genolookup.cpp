//
// C++ Implementation: genolookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genolookup.h"
#include <assert.h>

namespace Utility {

GenoLookup::~GenoLookup()
{
}
/**
 * Returns the mapped value setting it if necessary
 */
int GenoLookup::GetValue(const char* key)	{
	int value = 0;
	size_t size = strlen(key);
	char *lhs = new char[size+5];
	char *rhs = new char[size+5];
	const char *delPos = strchr(key, delim);
	
	if (delPos) {
		strncpy(lhs, key, delPos-key);
		lhs[delPos-key]='\0';
		strcpy(rhs, delPos+1);
	
		if (atoi(lhs) * atoi(rhs) == 0)
			value = NODATA;
		//Homozygous
		else if (strcmp(lhs, rhs) == 0)
			value = GetHomoZygote(key);
		else
			value = GetHeteroZygote(key);
	}
	else	{
		cout<<"Unbalanced Genotype: "<<key<<". Expected delimeter: "<<delim<<".\n";
		assert(0);
	}

	delete[] lhs;
	delete[] rhs;
	return value;
}

string GenoLookup::GetValue( int idx) {
	if (idx==NODATA)
		return "M";
	else if (idx==HOMOZYGOTE1) 
		return homozygote1;
	else if (idx==HETEROZYGOTE)
		return heterozygote1;
	else if (idx==HOMOZYGOTE2)
		return homozygote2;
	else
		return "M";
}
int GenoLookup::GetHeteroZygote(const char* value) {
	int returnVal = NODATA;

	if (strcmp(heterozygote1.c_str(), value) == 0)
		returnVal = HETEROZYGOTE;
	else if (strcmp(heterozygote2.c_str(), value) == 0)
		returnVal = HETEROZYGOTE;
	else if (heterozygote1 == "")	{
		returnVal = HETEROZYGOTE;
		heterozygote1 = value;
	} else if (heterozygote2 == "") {
		returnVal = HETEROZYGOTE;
		heterozygote2 = value;
	}
	else 	{
		nodata.push_back(value);
		printf("Skipping genotype: %s\n", value);
	}
	return returnVal;
		

}
int GenoLookup::GetHomoZygote(const char* value) {
	int returnValue = NODATA;

	if (strcmp(homozygote1.c_str(), value) == 0)
		returnValue = HOMOZYGOTE1;
	else if (strcmp(homozygote2.c_str(), value) == 0)
		returnValue = HOMOZYGOTE2;
	else if (homozygote1 == "")	{
		returnValue = HOMOZYGOTE1;
		homozygote1 = value;
	} else if (homozygote2 == "") {
		returnValue = HOMOZYGOTE2;
		homozygote2 = value;
	}
	else	{
		nodata.push_back(value);
		printf("Skipping genotype: %s\n", value);
	}
	return returnValue;
}
 ostream &GenoLookup::operator>>(ostream& os) {
	Report(&os);
	return os;
}
ostream *GenoLookup::Report(ostream* os)		{		///<Appends the mapping to the stream
	*os<<"\n***** "<<headerName<<" *****\n";
	
	*os<<"\t"<<homozygote1<<" -> Homozygote 1("<<HOMOZYGOTE1<<")\n";
	*os<<"\t"<<homozygote2<<" -> Homozygote 2("<<HOMOZYGOTE2<<")\n";
	*os<<"\t"<<heterozygote1<<", "<<heterozygote2<<" -> Heterozygotes("<<HETEROZYGOTE<<")\n";

	*os<<"\n";
	return os;
}	


}
