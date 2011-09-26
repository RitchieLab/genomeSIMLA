//
// C++ Interface: appconfiguration
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERAPPCONFIGURATION_H
#define BIOFILTERAPPCONFIGURATION_H

#include "utility/lineparser.h"

namespace Paris {

/**
@Brief Reads/writes Configuration and initializes application settings accordingly

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class AppConfiguration : public Utility::FileToMapExplicit {
public:
    AppConfiguration();
    ~AppConfiguration();

	/**
	 * @Brief Write the configuration to stream
	 */
	void WriteConfiguration(std::ostream& os);

	/**
	 * @Brief Initialize configuration with defaults
	 */
	void Init();

	/**
	 * @Brief Updates application data according to values from the configuration
	 */
	void ExecuteConfiguration();

	/**
	 * @Brief Produce a report header listing key settings from the configuration
	 */
	void ReportConfiguration(std::ostream& os);

	void PrintSet(const char *key, std::vector<std::string>& settings);
};



}

#endif
