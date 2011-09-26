/* 
 * File:   main.h
 * Author: torstees
 *
 * Created on January 6, 2010, 12:54 PM
 */

#ifndef _MAIN_H
#define	_MAIN_H

#include "parisapp.h"
#include "appconfiguration.h"
#include <map>
#include "utility/types.h"

namespace Paris {
class Main {
public:
	Main();
	Main(const Main& orig);
	virtual ~Main();

	bool ParseCmdLine(int argc, char **argv);
	int ParseCmd(int curr, int argc, char **argv);
	void PrintHelp();										///< Display usage details
	void PrintBanner();									///< Display details about the software
	AppConfiguration *LoadConfiguration(const char *cfgFilename);
	std::string GetReportPrefix();
	void RunCommands();
	void LoadDataPoints(std::vector<ParisApp::DataPoint>& data, std::multimap<std::string, uint >& snps);
	Utility::StringArray LoadGroupFile(const char *groups);
private:
	AppConfiguration cfg;								///< Configuration object
	std::string configFilename;						///< Configuration filename
	ParisApp parisApp;									///< This does all the work

	struct ParisAction {
		enum Action {
			NoAction,										///< No particular action. This is the default
			ParseError,										///< Problems were found
			PrintSampleConfig,							///< PrintSampleConfiguration
			InvestigatePathways,							///< Group Investigation Report
			ListGroups,										///< List available groups
			ListPopulationIDs,							///< List Population IDs
			ProcessReport,									///<
			WriteBinReport									///< Report Bin# and lowest pvalue for each feature
		};
	};
	ParisAction::Action action;
};
}
#endif	/* _MAIN_H */

