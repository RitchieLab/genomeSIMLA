//
// C++ Interface: configuration
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESECONFIGURATION_H
#define ESECONFIGURATION_H

#include "utility/utility.h"
#include <iostream>

namespace ESE {

using namespace std;
using namespace Utility;

/**
@brief Basic configuration management for MDR-PDT

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Configuration : public Utility::ConfigurationParser
{
public:

	/**
	 * Inherited functionality from ConfigurationParser
	 */
	/**
 	 * Dumps a Text Report summarizing the configuration options encountered
	 */
	void GenerateReport(ostream& os);	
	void GenerateReport();
	/**
	 * Setup the key value assignments
	 */
	bool SetValues(const char *key, const char *value);

	bool Validate();
	
	/**
	 * Perform any data mods required after loading is completed
	 */
	void PostLoad();

	BasicLog *GetModelStream();
	BasicLog *GetMendelStream();
	BasicLog *GetGeneralStream();

	string errorMsg;								///<Error associated with parsing and setting up configuration
	string dataFile;								///<File containing genetic data
	uint genotypicErrorResponse;					///<Determines how the program responds to genotype errors 
	uint genotypicErrorThreshold;					///<Genotype error threshold required before the system will purge the family					
	string projectName;								///<The base filename for each of the output files
	string mendelExt;								///<The extension for outputting the mendel error report
	string modelsExt;								///<The extension associated with the model output file
	string generalExt;

	/**
	 * Singleton stuff
	 */
	static Configuration *Instance();
	static void Release();
protected:
	BasicLog *modelLog;
	uint modelLogLevel;
	
	BasicLog *mendelLog;
	uint mendelLogLevel;

	BasicLog *generalLog;
	uint generalLogLevel;

    Configuration();

    ~Configuration();
	static uint _instanceCount;
	static Configuration *_instance;

class ConfigKeys {
	public:	
	static string DataSourceFile;					///<The key associated with setting the data source 
	static string ErrorThreshold;					///<The key associated Genotype error threshold	
	static string ErrorResponse;					///<The key associated with error response to genotype errors

	static string BaseReportName;					///<The base part of the filename

	static string ExtErrorReport;					///<The extension associated with the mendal error report
	static string ExtModelsReport;					///<The extension associated with the models output
	static string ExtGeneralReport;					///<The extension associated with the general output
	static string ExtPTestReport;					///<The extension where the PTest results are to be written

};
};

inline
BasicLog *Configuration::GetGeneralStream() {
	if (generalLog)
		return generalLog;
	else if (strlen(generalExt.c_str())==0 || strcmp(generalExt.c_str(), "STDOUT")==0) {
		generalLog = new AsciiLog("STDOUT", generalLogLevel);
	}
	else {
		string filename = projectName + "." + generalExt;
		generalLog = new AsciiLog(filename.c_str(), generalLogLevel);
	}
	return generalLog;
}

inline
BasicLog *Configuration::GetMendelStream() {
	if (mendelLog)
		return mendelLog;
	else if (strlen(mendelExt.c_str())==0 || strcmp(mendelExt.c_str(), "STDOUT")==0) {
		mendelLog = new AsciiLog("STDOUT", mendelLogLevel);
	}
	else {
		string filename = projectName + "." + mendelExt;
		mendelLog = new AsciiLog(filename.c_str(), mendelLogLevel);
		mendelLog->Write("Mendelian Error Log:\n");
	}
	return mendelLog;
}

inline
BasicLog *Configuration::GetModelStream() {
	if (modelLog)
		return modelLog;
	else if (strlen(modelsExt.c_str())==0 || strcmp(modelsExt.c_str(), "STDOUT")==0) {
		modelLog = new AsciiLog("STDOUT", modelLogLevel);
	}
	else {
		string filename = projectName + "." + modelsExt;
		modelLog = new AsciiLog(filename.c_str(), modelLogLevel);
	}
	return modelLog;
}

inline
Configuration::Configuration()
 : ConfigurationParser(), modelLog(NULL), modelLogLevel(10), mendelLog(NULL), mendelLogLevel(10), generalLog(NULL), generalLogLevel(10) {
}

inline
Configuration::~Configuration()
{
	if (generalLog)
		delete generalLog;
	if (modelLog)
		delete modelLog;
	if (mendelLog)
		delete mendelLog;
}

inline
Configuration *Configuration::Instance() {
	if (_instance == NULL)
		_instance = new Configuration();
	_instanceCount++;
	return _instance;
}

inline
void Configuration::Release() {
	if (--_instanceCount < 1) {
		delete _instance;
		_instance = NULL;
	}
}

}

#endif
