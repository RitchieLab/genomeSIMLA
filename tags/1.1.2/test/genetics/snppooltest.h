//
// C++ Interface: snppooltest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_TESTSNPPOOLTEST_H
#define GENETICS_TESTSNPPOOLTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "genetics/snppool.h"


namespace Genetics {

namespace Test {

/**
Test for SnpPool objects

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpPoolTest: public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( SnpPoolTest );
	CPPUNIT_TEST( TestAcquisitionAndRelease );
	CPPUNIT_TEST_SUITE_END();
public:
    SnpPoolTest();

    ~SnpPoolTest();

	void setUp();
	void tearDown();

	void TestAcquisitionAndRelease();
	
protected:
	SnpPool *pool;

};

}

}

#endif
