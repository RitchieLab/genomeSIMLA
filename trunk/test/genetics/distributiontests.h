//
// C++ Interface: distributiontests
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_TESTDISTRIBUTIONTESTS_H
#define GENETICS_TESTDISTRIBUTIONTESTS_H

#include <cppunit/extensions/HelperMacros.h>

namespace Genetics {

namespace Test {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class DistributionTests : public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( DistributionTests );
	CPPUNIT_TEST( TestNTestDistribution );
	CPPUNIT_TEST_SUITE_END();
public:
    DistributionTests();
    ~DistributionTests();

	void TestNTestDistribution();

};

}

}

#endif
