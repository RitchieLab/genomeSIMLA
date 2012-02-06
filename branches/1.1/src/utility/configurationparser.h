//
// C++ Interface: configurationparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYCONFIGURATIONPARSER_H
#define UTILITYCONFIGURATIONPARSER_H
#include "lineparser.h"
#include "types.h"
#include <boost/tokenizer.hpp>



namespace Utility {

using namespace boost;
//typedef tokenizer<char_separator<char> > cfgtokenizertype;



/**
@brief Baseclass for configuration parsers. 
Overload this class defining the GenerateReport and SetValues functions 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ConfigurationParser : public AsciiParser {
public:
    ConfigurationParser() {	}

    virtual ~ConfigurationParser() {}

	/**
 	 * Dumps a Text Report summarizing the configuration options encountered
	 */
	virtual void GenerateReport(ostream& os) = 0;	
	
	/**
	 * @brief Called to assign <i>value</i> to the entry, <i>key</i>
	 * @param key Left hand portion of the configuration line
	 * @param value Right hand portion of the configuration line
	 * @return true/false to indicate that it was a valid entry that was handled
	 */
	virtual bool SetValues(const char *key, const char *value, const char *remainder) = 0;

	/**
	 * @brief Check for missing data or impossible configuration settings
	 */
	virtual bool Validate() = 0;

	
	/**
	 * @brief Accept a line of input from an ascii text file
	 * @param line A single line of ascii text
	 * @return True if the line consisted of a valid entry
	 */
	bool ParseLine(const char *line, uint);

	/**
	 * @brief Derived classes will use this to extract boolean values from configuration files
	 */
	bool GetBoolean(const char *value);

	uint CountColumns(const char *line);
};


}

#endif
