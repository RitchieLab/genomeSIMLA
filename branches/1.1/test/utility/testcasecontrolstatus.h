//
// C++ Interface: testcasecontrolstatus
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITY_TESTTESTCASECONTROLSTATUS_H
#define UTILITY_TESTTESTCASECONTROLSTATUS_H

#include <cppunit/extensions/HelperMacros.h>
#include "utility/utility.h"
namespace Utility {

namespace Test {

/**
@brief tests associated with the CaseControlStatus object

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

using namespace std;


class TestCaseControlStatus : public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( TestCaseControlStatus );
	CPPUNIT_TEST( TestBalanced );
	CPPUNIT_TEST( TestImperfect );
	CPPUNIT_TEST( TestCopyConstructor );
	CPPUNIT_TEST( TestAppendStatus );
	CPPUNIT_TEST( TestSetStatus );
	CPPUNIT_TEST_SUITE_END();

public:
	TestCaseControlStatus() {}
	~TestCaseControlStatus() {}
	void TestBalanced();			///<Tests a set with balanced and complete data 
	void TestImperfect();			///<Tests a set with unbalanced and incomplete data
	void TestCopyConstructor();		///<Tests the assignment and propagation of data through copy contronstructor
	void TestSetStatus();
	void TestAppendStatus();
};


}

}

#endif
