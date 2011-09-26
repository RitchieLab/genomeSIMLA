//
// C++ Interface: snprepostest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESE_TESTSNPREPOSTEST_H
#define ESE_TESTSNPREPOSTEST_H

#include <cppunit/extensions/HelperMacros.h>
#include "genetics/snppool.h"
#include "genetics/snpaligned.h"

namespace ESE {

namespace Test {

using namespace Genetics;

/**
Test for the snp respository used in ESE

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpReposTest : public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( SnpReposTest );
	CPPUNIT_TEST( TestLoading );
	CPPUNIT_TEST( TestMDRLoading );
//	CPPUNIT_TEST( TestBinLoading );
	CPPUNIT_TEST( TestMultiLoading );
	CPPUNIT_TEST( TestPerfect );
	CPPUNIT_TEST( TestMaxDiff );
	CPPUNIT_TEST( TestSelection );
	CPPUNIT_TEST_SUITE_END();
	SnpPool *pool;
public:
    SnpReposTest();

    ~SnpReposTest();

	void TestLoading();
	void TestMultiLoading();
	void TestBinLoading();
	void TestPerfect();
	void TestMaxDiff();
	void TestSelection();
	void TestMDRLoading();
	void CompareSnps(SnpAligned *lhs, SnpAligned *rhs);


};

}

}

#endif
