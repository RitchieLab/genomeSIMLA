//
// C++ Interface: lineparsing
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITY_TESTLINEPARSING_H
#define UTILITY_TESTLINEPARSING_H

#include <cppunit/extensions/HelperMacros.h>

namespace Utility {

namespace Test {

/**
@brief Tests functionality of the line parsing class as well as any related functions/classes

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LineParsing: public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( LineParsing );
	CPPUNIT_TEST( TestExclusionList );
	CPPUNIT_TEST_SUITE_END();
public:
    LineParsing();

    ~LineParsing();

	void TestExclusionList();

};

}

}

#endif
