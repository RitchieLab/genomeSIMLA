//
// C++ Implementation: application
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "application.h"
#include <iostream>
#include <fstream>

namespace Utility {

Application::Application() 
{
	major=0;
	minor=1;
	bugFixes=0;
	cmdError = false;
}


Application::~Application()
{
}

void Application::PrintBanner() {
	printf("%s - %s: %d.%d.%d build %d (%s)\n", appname.c_str(), appfunction.c_str(), major, minor, bugFixes, buildNumber, buildDate.c_str());
	printf("Vanderbilt University\nCenter for Human Genetics Research\n");
	printf("%s\n\n", authors.c_str());
}


void Application::PrintBanner(BasicLog *log) {
	log->Write(0, "%s - %s: %d.%d.%d build %d (%s)\n", appname.c_str(), appfunction.c_str(), major, minor, bugFixes, buildNumber, buildDate.c_str());
	log->Write( 0, "Vanderbilt University\nCenter for Human Genetics Research\n");
	log->Write( 0, "%s\n\n", authors.c_str());
}


uint Application::HandleArgument(char *flag, char *value)	{
	printf("Unknown command: %s %s\n", flag, value);
	cmdError = true;
	return 1;
} 


bool Application::ParseCmdLine(int argc, char **argv) {
	int current = 1;

	while (!cmdError && current < argc) 	{
		if (argc >= current+1)
			current += HandleArgument(argv[current], argv[current+1]);
	}

	return !cmdError && argc > 1;
}





}


#ifdef TEST_APP
#include <gtest/gtest.h>
#include "utility/random.h"



int
main( int argc, char* argv[] )
{
	Utility::Random::globalGenerator = Utility::Random(1371);
	//Utility::Random rnd = Utility::Random(1323);
::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}

#endif
