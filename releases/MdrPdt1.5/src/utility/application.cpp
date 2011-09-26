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
	printf("%s - %s: %d.%d.%d \n", appname.c_str(), appfunction.c_str(), major, minor, bugFixes);
	printf("Vanderbilt University\nCenter for Human Genetics Research\n");
	printf("%s\n\n", authors.c_str());
}

void Application::PrintBanner(BasicLog *log) {
	log->Write( 0, "%s - %s: %d.%d.%d\n", appname.c_str(), appfunction.c_str(), major, minor, bugFixes);
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
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>


int
main( int argc, char* argv[] )
{
  // Create the event manager and test controller
  CPPUNIT_NS::TestResult controller;

  // Add a listener that colllects test result
  CPPUNIT_NS::TestResultCollector result;
  controller.addListener( &result );        

  // Add a listener that print dots as test run.
  CPPUNIT_NS::BriefTestProgressListener progress;
  controller.addListener( &progress );      

  // Add the top suite to the test runner
  CPPUNIT_NS::TestRunner runner;
  runner.addTest( CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest() );
  runner.run( controller );

  // Print test in a compiler compatible format.
  CPPUNIT_NS::CompilerOutputter outputter( &result, std::cerr );
  outputter.write(); 

  return result.wasSuccessful() ? 0 : 1;
}

#endif
