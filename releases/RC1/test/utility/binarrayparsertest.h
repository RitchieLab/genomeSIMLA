//
// C++ Interface: binarrayparsertest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITY_TESTBINARRAYPARSERTEST_H
#define UTILITY_TESTBINARRAYPARSERTEST_H

#include <cppunit/extensions/HelperMacros.h>

namespace Utility {

namespace Test {

/**
Test for BinArray Functionality

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BinArrayParserTest: public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( BinArrayParserTest );
	CPPUNIT_TEST( LearnBitManipulation );
	CPPUNIT_TEST( DoTests );
	CPPUNIT_TEST_SUITE_END();

public:
    BinArrayParserTest();

    ~BinArrayParserTest();

	void LearnBitManipulation();
	void DoTests();


	unsigned char GetBinWindow(unsigned char val, int frame, int framesize, int length, int offset=0);

};

}

}

#endif
