//
// C++ Interface: application
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef APPLICATION_H
#define APPLICATION_H

#include <string>
#include "types.h"
#include "basiclog.h"

namespace Utility {

	using namespace std;

/**
@brief Base class for various applications. This provides some simple functions and a generic interface to follow.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Application{
public:
    Application();

    virtual ~Application();

	/**
 	 * @brief This will let the program set up the necessary parameters based on cmd line arguments. 
	 * @return true indicates that execution can begin
 	 */
	virtual bool ParseCmdLine(int argc, char** argv);

	/**
	 * @brief Prints the help contents
	 */
	virtual void PrintHelp() = 0;
	
	/**	
	 * @brief Starts execution
	 */
	virtual void Start() = 0;

	/**
	 * @brief Give the app a chance to run a quick sanity check
	 */
	virtual bool PerformSanityCheck() { return true; }

	/**
	 * @brief Displays basic information about the application
	 */
	void PrintBanner();
	void PrintBanner(BasicLog *log);
protected:
	/**
	 * @brief Parses individual commands and executes necessary tasks associated with them
	 */
	virtual uint HandleArgument(char *flag, char *value);

	/**
	 * @brief Do a quick sanity check on the settings provided
	 */
	bool VerifyConfiguration();

	bool cmdError;						///<Indicates that an error was found with the command line arguments

	string appname;						///<This is the program's name. Simply used as part of the generic banner and help stuff
	string appfunction;					///<Quick overview of the application's purpose
	string authors;	

	uint major;							///<Some versioning control information
	uint minor;							///<Some versioning control information
	uint bugFixes;						///<Some versioning control information
	uint buildNumber;					
	string buildType;					///<Debug, if it's debug mode
	string buildDate;					///<Date of build


};

}

#endif
