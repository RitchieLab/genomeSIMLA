//
// C++ Interface: stringlookup
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILSTRINGLOOKUP_H
#define UTILSTRINGLOOKUP_H

#include <map>
#include <iostream>
#include <string>
#include "utility.h"

namespace Utility {

using namespace std;

typedef map<string, uint> LookupType;
/**
@brief Quick lookup table that allows the user to dump the mapping between the strings and the values to a stream.


	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class StringLookup{
public:
    StringLookup();				
    ~StringLookup();
	/**
	 * @brief Sets the name as it should appear in the report
	 * @param name The name of the lookup table
	 */
	void SetHeader(const char* name);		
	/**
	 * @brief Returns the mapped value setting it if necessary
	 * @param key The value that is being observed
	 * If this value doesn't exist, it will create a new index to it. Otherwise, it returns the 
	 * previously observed index.
	 */
	uint GetValue(const char* key);			
	/**
	 * @brief Appends the mapping to the stream	
	 * @param os The stream to be written to
	 */
	void Report(ostream& os);				
	/**
	 * @brief Returns the number of entries in the table
	 */
	uint Size();							
protected:
	LookupType lookupTable;					///<This is the lookup table itself
	uint nextVal;							///<This keeps up with the values as they come in
	string headerName;						///<This is the name as it appears in the report

};



}

#endif
