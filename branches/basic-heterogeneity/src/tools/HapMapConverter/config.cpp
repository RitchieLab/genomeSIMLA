//
// C++ Implementation: config
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <iostream>
#include <iomanip>
#include "config.h"
#include "utility/utility.h"
#include <assert.h>
#include <sstream>
#include "utility/strings.h"
#include <string>
#include <map>

namespace Tools {

const char *Config::ConfigKeywords::SampleInputFile="SAMPLE_INPUTFILE";
const char *Config::ConfigKeywords::PhaseInputFile="PHASE_INPUTFILE";
const char *Config::ConfigKeywords::LegendInputFile="LEGEND_INPUTFILE";
const char *Config::ConfigKeywords::InclusionList="INCLUSION_INPUTFILE";


using namespace std;

Config::Config() : configurationFilename(""), inclusionsList("")  {	
	destinationFilename = configurationFilename + string("dataset.ped");
}
Config::~Config()	{ }

void Config::SetConfig(const char *cfgFilename) {
	configurationFilename = cfgFilename;
	destinationFilename = configurationFilename + string("dataset.ped");
}

bool Config::SetValues(const char *key, const char *value, const char *line) {
	bool success=true;
	assert(strlen(key)>0);


	if (strcmp(key, ConfigKeywords::SampleInputFile) == 0)  
		sampleInputFile = value;
	else if (strcmp(key, ConfigKeywords::PhaseInputFile) == 0)
		phaseInputFile = value;
	else if (strcmp(key, ConfigKeywords::LegendInputFile) == 0)
		legendInputFile = value;
	else if (strcmp(key, ConfigKeywords::InclusionList) == 0)
		inclusionsList = value;
	else {
		cout<<"!! Unhandled line: "<<line<<"\n";
		success = false;
	}
	return success;
}


string Config::GetRemainderString(const char *line) {
	string junk;
	string var;
	
	stringstream ss(line);
	ss>>junk;

	while (!ss.eof())  {
		ss>>junk;
		var += " " + junk;
		junk = "";
	}
	
	return var;

}

void Config::GenerateReport( ostream &os) {
	uint width = 45;
	cout<<"Configuration Settings:\n";
	cout<<setw(width)<<"Sample Source: "<<sampleInputFile<<"\n";
	cout<<setw(width)<<"Phased Source: "<<phaseInputFile<<"\n";
	cout<<setw(width)<<"Legend Source: "<<legendInputFile<<"\n";
	cout<<setw(width)<<"Inclusion List: "<<inclusionsList<<"\n";

}

/**
 * @brief Check for missing data or impossible configuration settings
 */
bool Config::Validate() {
	bool isValid = true;

	

	return isValid;
}




}
