//
// C++ Interface: fileparsertests
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_TESTFILEPARSERTESTS_H
#define GENETICS_TESTFILEPARSERTESTS_H

#include <cppunit/extensions/HelperMacros.h>
#include "genetics/snppool.h"
#include "genetics/snpaligned.h"

namespace Genetics {

namespace Test {

/**
@brief test for different file parsers

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class FileParserTests: public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( FileParserTests );
	CPPUNIT_TEST( TestMdrLineParser );
	CPPUNIT_TEST( TestMdrLoading);
	CPPUNIT_TEST_SUITE_END();
public:
    FileParserTests();
    ~FileParserTests();

	void CompareSnps(SnpAligned *lhs, SnpAligned *rhs);

	void TestMdrLineParser();
	void TestMdrLoading();

protected:
	SnpPool *pool;
};


}

}

#endif
