//
// C++ Interface: config
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef CONFIG_H
#define CONFIG_H

#include "utility/utility.h"

namespace Tools {

using namespace std;
/**
@brief The configuration parser for genomeSIM

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Config : public Utility::ConfigurationParser {
public:
    Config();

    ~Config();

	/**
 	 * Dumps a Text Report summarizing the configuration options encountered
	 */
	virtual void GenerateReport(ostream& os);	
	
	/**
	 * @brief Called to assign <i>value</i> to the entry, <i>key</i>
	 * @param key Left hand portion of the configuration line
	 * @param value Right hand portion of the configuration line
	 * @return true/false to indicate that it was a valid entry that was handled
	 */
	virtual bool SetValues(const char *key, const char *value, const char *remainder);

	/**
	 * @brief Check for missing data or impossible configuration settings
	 */
	virtual bool Validate();

	string GetRemainderString(const char *line);

	void SetConfig(const char *cfgFilename);

	struct ConfigKeywords {
		static const char *SampleInputFile;
		static const char *PhaseInputFile;
		static const char *LegendInputFile;
		static const char *InclusionList;
	
	};

	/**
	 * @brief Source filename
	 */
	string configurationFilename;

	/**
	 * @brief the ids assocaited with the individuals
	 */
	string sampleInputFile;

	/**
	 * @brief The phased data
	 */
	string phaseInputFile;

	/**
	 * @brief the rs-number information
	 */
	string legendInputFile;

	/**
	 * @brief the name of the file we create at the very end
	 */
	string destinationFilename;

	/**
	 * @brief list of snps to be included in the output file
	 */
	string inclusionsList;
};

					

}

#endif
