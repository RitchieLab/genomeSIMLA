//
// C++ Interface: genolookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GAWGENOLOOKUP_H
#define GAWGENOLOOKUP_H

#include "typelookup.h"
#include "types.h"
#include <string>
#include <vector>
#include <assert.h>

/*
#define HOMOZYGOTE1 	0
#define HETEROZYGOTE 	1
#define HOMOZYGOTE2 	2
#define NODATA      	-1
*/
#define HOMOZYGOTE1 	1
#define HETEROZYGOTE 	2
#define HOMOZYGOTE2 	3
#define NODATA      	0


namespace Utility {


using namespace std;

/**
A specialized lookup table used to correctly set up and retain the mapping details for allelic coded genotypes.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenoLookup : public TypeLookup
{
public:
    GenoLookup();
    ~GenoLookup();

	void SetHeader(const char* name);		///<Sets the name as it should appear in the report
	string GetValue( int idx);
	virtual int GetValue(const char* key);	///<Returns the mapped value setting it if necessary
	ostream *Report(ostream* os);			///<Appends the mapping to the stream
	const char *GetHeader();				///<Returns the header string
	ostream &operator>>(ostream& os);		///<Feed the contents of this structure to the stream
	void SetDelimeter(char del);			///<Sets the delimeter used to parse values from a file
	GenoLookup &operator=(const GenoLookup&o);		

	/**
	 * @brief Parses the binary representation of the genotype data
	 * @param encodedGT This is the encoded heterozygote
	 * @note Based on Lance's genotype encoding: (it's an 8 bit unsigned char) xxAAABBB
	 * <UL><LI>000 - M
	 * <LI>001 - A
	 * <LI>010 - C
 	 * <LI>011 - G
	 * <LI>100 - T
	 * </UL>
	 */
	bool ParseBinaryGenotypes(unsigned char encodedGT);
	unsigned char GetEncodedGenotypes();

	const char*GetHeterozygote();
	const char*GetHomozygote1();
	const char*GetHomozygote2();

	void Reset();
protected:
	virtual int GetHeteroZygote(const char *key);
	virtual int GetHomoZygote(const char *key);

	char delim;								///<This is the string containing the delimeter used to separate the two alleles
	string headerName;						///<This is the name as it appears in the report

	string homozygote1;						///<The first homozygote encountered
	string homozygote2;						///<The second homozygote encountered
	string heterozygote1;					///<One of the heterozygotes
	string heterozygote2;					///<The other one
	char SetupHomoZygotes(string& homoz, uint coding);
	void SetupHeteroZygotes(char hz1, char hz2);
	vector<string> nodata;					///<To keep up with the entries that were found with no data
};

inline
void GenoLookup::Reset() {
	homozygote1=homozygote2=heterozygote1=heterozygote2="";
	headerName="";
}
inline
const char *GenoLookup::GetHeterozygote() {
	const char *htz=heterozygote1.c_str();
	if (strcmp(htz, "")==0)
		heterozygote1="Aa";
	return heterozygote1.c_str();
}

inline
const char *GenoLookup::GetHomozygote1() {
	const char *hmz=homozygote1.c_str();
	if (strcmp(hmz, "")==0)
		homozygote1="AA";
	return homozygote1.c_str();
}

inline
const char *GenoLookup::GetHomozygote2() {
	const char *hmz=homozygote2.c_str();
	if (strcmp(hmz, "")==0)
		homozygote2="aa";
	return homozygote2.c_str();
}

inline
bool GenoLookup::ParseBinaryGenotypes(unsigned char encodedGT) {
	char hz2 = SetupHomoZygotes(homozygote2, encodedGT&7);
	char hz1 = SetupHomoZygotes(homozygote1, encodedGT>>3);
	SetupHeteroZygotes(hz1, hz2);
	//Set this up to be unknown
//	else
//		heterozygote1 = heterozygote2 = homozygote2 = homozygote1="MM";
	return true;
}

inline
void GenoLookup::SetupHeteroZygotes(char hz1, char hz2) {
	char hz[3];
	hz[0]=hz1;
	hz[1]=hz2;
	hz[2]='\0';
	heterozygote1=hz;
	hz[0]=hz2;
	hz[1]=hz1;
	heterozygote2=hz;
}	
inline
char GenoLookup::SetupHomoZygotes(string& homoz, uint coding) {
	static char alleles[6]="MACGT";
	char hz[3];
	assert(coding>0 && coding<5);
	char allelicValue=alleles[coding];
	hz[0]=allelicValue;
	hz[1]=allelicValue;
	hz[2]='\0';
	homoz=hz;
	return allelicValue;
}
	

inline
GenoLookup &GenoLookup::operator=(const GenoLookup &o) {
	delim=o.delim;
	headerName=o.headerName;
	homozygote1=o.homozygote1;
	homozygote2=o.homozygote2;	
	heterozygote1=o.heterozygote1;
	heterozygote2=o.heterozygote2;
	return *this;
}


inline
GenoLookup::GenoLookup() : delim('/'), headerName(""), homozygote1(""), homozygote2(""), heterozygote1(""), heterozygote2("") {}


inline
void GenoLookup::SetDelimeter(const char del) {
	delim=del;
}

inline
void GenoLookup::SetHeader(const char *name) {
	headerName = name;
}

inline
const char *GenoLookup::GetHeader() {
	return headerName.c_str();
}




}

#endif
