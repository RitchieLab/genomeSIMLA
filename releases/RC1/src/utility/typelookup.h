//
// C++ Interface: typelookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYTYPELOOKUP_H
#define UTILITYTYPELOOKUP_H

#include <iostream>

namespace Utility {

using namespace std;

/**
Base class for the various genotype and phenotype lookups

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class TypeLookup{
public:
    TypeLookup();

    virtual ~TypeLookup();


	void SetHeader(const char* name);			///<Sets the name as it should appear in the report
	virtual int GetValue(const char* key)	=0;	///<Returns the mapped value setting it if necessary
 	virtual string GetValue(int i)			=0;
	virtual ostream *Report(ostream* os)	=0;	///<Appends the mapping to the stream
	virtual const char *GetHeader()			=0;
	virtual ostream &operator>>(ostream& os)=0;	
	
	int GetNotEncodedIdx();
	void SetNotEncodedIdx(int idx);
protected:
	string headerName;							///<This is the name as it appears in the report
	int badIdx;

};

inline
int TypeLookup::GetNotEncodedIdx() {
	return badIdx;
}

inline
void TypeLookup::SetNotEncodedIdx( int idx) {
	badIdx=idx;
}

inline
void TypeLookup::SetHeader(const char *name) {
	headerName = name;
}

inline
TypeLookup::TypeLookup() : headerName(""), badIdx(-1) {}

inline
TypeLookup::~TypeLookup() {}

}

#endif
