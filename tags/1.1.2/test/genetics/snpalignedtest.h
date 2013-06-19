//
// C++ Interface: snpalignedtest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TESTSNPALIGNEDTEST_H
#define TESTSNPALIGNEDTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "genetics/snppool.h"


namespace Genetics {
namespace Test {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpAlignedTest : public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( SnpAlignedTest );
	CPPUNIT_TEST( TestBitSetAssignment);
	CPPUNIT_TEST( TestLabels);
	CPPUNIT_TEST( TestUtil );
	CPPUNIT_TEST( TestImporting );
	CPPUNIT_TEST( ImportCompressedBinary );
	CPPUNIT_TEST( TestPunett );
	CPPUNIT_TEST( TestMerge );
	CPPUNIT_TEST( TestGenoLookup );
	CPPUNIT_TEST( TestDescriptors );
	CPPUNIT_TEST_SUITE_END();
public:
    SnpAlignedTest();

    ~SnpAlignedTest();
	
	/**
	 * Called once during testing to set up and tear down structures
	 */
	void setUp();
	void tearDown();

	void TestUtil();
	void TestImporting();
	void TestPunett();
	void TestMerge();
	void ImportCompressedBinary();
	void TestLabels();
	void TestDescriptors();
	void TestGenoLookup();

	void TestBitSetAssignment();
	void TestAssignment();

protected:
	SnpPool *pool;
};

}
}

#endif
