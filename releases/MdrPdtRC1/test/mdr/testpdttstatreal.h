//
// C++ Interface: testpdttstatreal
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESE_TESTTESTPDTTSTATREAL_H
#define ESE_TESTTESTPDTTSTATREAL_H
#include "testpdttstat.h"

namespace MDRPDT {

namespace Test {

/**
@brief Perform tests on real data and verify their outcome. These tests shouldn't ever be delivered outside of the lab, since the data isn't necessarily data we have the right to distribute.

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class TestPdtTStatReal : public TestPdtTStat
{
	CPPUNIT_TEST_SUITE( TestPdtTStatReal );
	CPPUNIT_TEST( TestModelsALZ );
	CPPUNIT_TEST( TestModelsSCZ );
	CPPUNIT_TEST( TestModelsBPADCIT );
	CPPUNIT_TEST_SUITE_END();
public:
    TestPdtTStatReal();
    ~TestPdtTStatReal();

	void TestModelsSCZ();
	void TestModelsALZ();
	void TestModelsBPADCIT();

	void TestModelAlz();

};

}

}

#endif
