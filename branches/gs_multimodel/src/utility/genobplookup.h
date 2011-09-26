//
// C++ Interface: genobplookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYGENOBPLOOKUP_H
#define UTILITYGENOBPLOOKUP_H
#include "genolookup.h"
namespace Utility {

/**
Provides the ability to parse lines of Base Pair encoded genetypes

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenoBPLookup : public GenoLookup
{
public:
    GenoBPLookup();

    ~GenoBPLookup();

	int GetValue(const char* key);			///<Returns the mapped value setting it if necessary
protected:
	int GetHomoZygote(const char* value);

};

inline
GenoBPLookup::GenoBPLookup()
{
}

inline
GenoBPLookup::~GenoBPLookup()
{
}



}

#endif
