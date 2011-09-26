//
// C++ Implementation: configuration
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "configuration.h"


#define MDRPDT_THRESHOLD_LEVEL_MAX 3

namespace MDRPDT {
uint Configuration::_instanceCount 		= 0;
Configuration *Configuration::_instance	= NULL;
string Configuration::ConfigKeys::DataSourceFile 		= "InputFile";			///<The key associated with setting the data source 
string Configuration::ConfigKeys::ErrorThreshold		= "Threshold";			///<The key associated Genotype error threshold	
string Configuration::ConfigKeys::ErrorResponse			= "ErrorResponse";		///<The key associated with error response to genotype errors
string Configuration::ConfigKeys::BaseReportName		= "BaseReportName";		///<The key associated with setting the repot names
string Configuration::ConfigKeys::ExtErrorReport 		= "ExtensionMendelErrors";	///<The extension associated with the mendal error report
string Configuration::ConfigKeys::ExtModelsReport		= "ExtensionModelsReport";	///<The extension associated with the models output
string Configuration::ConfigKeys::ExtGeneralReport		= "ExtensionGeneralOutput";	///<The extension associated with the general output


/**
 * Dumps a Text Report summarizing the configuration options encountered
 */
void Configuration::GenerateReport(ostream& os) {
	string response[] = {"Report Only", "Strip Invalid from affected children only", "Strip invalid loci from all children in family", "Remove family if error count for a single child exceeds threshold else strip ivalid loci from all"};
	os<<"Inputfile               : "<<dataFile<<"\n";
	os<<"Genotype Error Response : \n";

	for (uint i=0; i<4; i++) {
		if (genotypicErrorResponse == i)
			os<<"\t* "<<response[i]<<"\n";
		else
			os<<"\t  "<<response[i]<<"\n";
	}
	if (genotypicErrorResponse > 3) 
		os<<"\t* Invalid Selection\n";

	os<<"Genotype Error Threshold : "<<genotypicErrorThreshold<<"\n";
}

void Configuration::GenerateReport() {
	ostream *os = GetGeneralStream()->GetStream();


	string response[] = {"Report Only", "Strip Invalid from affected children only", "Strip invalid loci from all children in family", "Remove family if error count for a single child exceeds threshold else strip ivalid loci from all"};
	*os<<"Inputfile               : "<<dataFile<<"\n";
	*os<<"Genotype Error Response : \n";

	for (uint i=0; i<4; i++) {
		if (genotypicErrorResponse == i)
			*os<<"\t* "<<response[i]<<"\n";
		else
			*os<<"\t  "<<response[i]<<"\n";
	}
	if (genotypicErrorResponse > 3) 
		*os<<"\t* Invalid Selection\n";

	*os<<"Genotype Error Threshold : "<<genotypicErrorThreshold<<"\n";
	
	*os<<"Medelian Errors:           "<<GetMendelStream()->GetFilename()<<"\n";
	*os<<"Model Report:              "<<GetModelStream()->GetFilename()<<"\n";
	*os<<"General Output:            "<<GetGeneralStream()->GetFilename()<<"\n";

}

/**
 * Setup the key value assignments
 */
bool Configuration::SetValues(const char *key, const char *value) {
	bool success=true;
	assert(strlen(key)>0);
	if (strcmp(key, ConfigKeys::DataSourceFile.c_str()) == 0) 
		dataFile = value;
	else if (strcmp(key, ConfigKeys::ErrorThreshold.c_str()) == 0)
		genotypicErrorThreshold=atoi(value);
	else if (strcmp(key, ConfigKeys::ErrorResponse.c_str()) == 0)
		genotypicErrorResponse=atoi(value);
	else if (strcmp(key, ConfigKeys::BaseReportName.c_str()) == 0)
		projectName=value;
	else if (strcmp(key, ConfigKeys::ExtErrorReport.c_str()) == 0)
		mendelExt=value;
	else if (strcmp(key, ConfigKeys::ExtModelsReport.c_str()) == 0)
		modelsExt=value;
	else if (strcmp(key, ConfigKeys::ExtGeneralReport.c_str()) == 0)
		generalExt=value;
	return success;
}

bool Configuration::Validate() {
	bool success = false;
	if (genotypicErrorResponse > MDRPDT_THRESHOLD_LEVEL_MAX)
		errorMsg="Threshold invalid";
	

	else 
		success= true;
	return success;
}


void Configuration::PostLoad() {

}

}
