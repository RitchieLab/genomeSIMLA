//
// C++ Interface: AppConfig
//
// Description: 
// Configuration object for MdrPDT (version 2) 
// This class exposes a single global object applicationCfg and 
// has the ability to read/write configurations 
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef APP_CONFIG_H
#define APP_CONFIG_H

#include "utility/lineparser.h"
#include <iostream>
#include <string>

namespace MdrPDT {

class ApplicationConfiguration;


class ApplicationConfiguration : public Utility::FileToMapExplicit {
public:
	virtual void WriteConfiguration(std::ostream& os)=0;
	virtual void Init()=0;
	virtual void ExecuteConfiguration()=0;
	virtual void ReportConfiguration(std::ostream& os)=0;
};



class PedigreeConfiguration : public ApplicationConfiguration {
	/**
	 * @brief Initialize with default values
	 */
	void Init();

	/**
	 * @brief Write out the contents to the file
	 */
	void WriteConfiguration(std::ostream& os);
	void ExecuteConfiguration();

	void ReportConfiguration(std::ostream& os);
};



class ConfigurationReader {
public:
	ConfigurationReader() { }
	~ConfigurationReader();

	/**
	 * @brief Load the configuration file
	 */
	ApplicationConfiguration* Load(const char *filename = NULL);
	/**
	 * @brief Report the configuration details to os
	 */
	void ReportConfiguration(std::ostream& os) { appConfiguration->ReportConfiguration(os); }

	/**
	 * @brief This variable is open for all to use (created and destroyed by main)
	 */
	static ApplicationConfiguration *appConfiguration;

	/**
	 * @brief Retuns the filename associated with the current configuration
	 */
	std::string GetFilename() { return filename;}
protected:

	std::string filename;
};


}

#endif

