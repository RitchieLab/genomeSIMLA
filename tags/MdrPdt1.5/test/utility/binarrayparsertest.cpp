//
// C++ Implementation: binarrayparsertest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "binarrayparsertest.h"
#include "utility/binarrayparser.h"
#include "utility/utility.h"
#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION( Utility::Test::BinArrayParserTest );

using namespace std;

namespace Utility {

namespace Test {

BinArrayParserTest::BinArrayParserTest()
{
}


BinArrayParserTest::~BinArrayParserTest()
{
}

unsigned char BinArrayParserTest::GetBinWindow(unsigned char val, int frame, int framesize, int readlength, int offset/*=0*/) {
	int start=frame*framesize+offset;
	int rFallout=start+readlength;
				
	unsigned char v=val;	
	v=v<<start;
	v>>=start+8-rFallout;
	v<<=framesize-readlength;
	return v;
}
	
	

void BinArrayParserTest::LearnBitManipulation() {
	unsigned char test = 69;				//01000101 -> 010 001 010
	CPPUNIT_ASSERT(GetBinWindow(test, 0, 3, 3) == 2);
	CPPUNIT_ASSERT(GetBinWindow(test, 1, 3, 3) == 1);

	//If this were under normal circumstances, we would need to grab the first bit from the next integer
	unsigned char val =GetBinWindow(test, 2, 3, 2);
	CPPUNIT_ASSERT(val == 2);
	test = 174;								//10101110 -> 101 011 10
	val|=GetBinWindow(test, 0, 1, 1);
	CPPUNIT_ASSERT(val == 3);


	CPPUNIT_ASSERT(GetBinWindow(test, 0, 3, 3) == 5);
	CPPUNIT_ASSERT(GetBinWindow(test, 1, 3, 3) == 3);
	CPPUNIT_ASSERT(GetBinWindow(test, 2, 3, 2) == 4);

	//Lets try offsets
	CPPUNIT_ASSERT(GetBinWindow(test, 0, 3, 3, 1) == 2);
	CPPUNIT_ASSERT(GetBinWindow(test, 1, 3, 3, 1) == 7);
	CPPUNIT_ASSERT(GetBinWindow(test, 2, 3, 2, 1) == 0);	
}

void BinArrayParserTest::DoTests() {

	//Here is the data:
	// 2,  1, 2,   4,  5,  6,    1,  0,  2,  2,2,1,2,4,5,6,1,0,2,2
	//010,001,01 0,100,101,1 10,001,000, 010,010,010,001,010,100,101,110,001,000,010,010, 000,0(padding of 1 set and a partial)
	//01000101   01001011    10001000    01001001 00010101 00101110 00100001 00100000
	//69,          75           136         73       21       46      33         32
	unsigned char testdata[] = { 69,75,136,73,21,46,33,32};
	BinArrayParser<unsigned char> binParser(testdata, 8, 3);	

	CPPUNIT_ASSERT(binParser.GetCount() == 8);
	CPPUNIT_ASSERT(binParser.GetFrameSize() == 3);
	CPPUNIT_ASSERT(binParser.GetWordSize() == 8);
	
	//int idx = 0;

	CPPUNIT_ASSERT(binParser.GetNextValue()==2);
	CPPUNIT_ASSERT(binParser.GetNextValue()==1);
	CPPUNIT_ASSERT(binParser.GetNextValue()==2);
	CPPUNIT_ASSERT(binParser.GetNextValue()==4);
	CPPUNIT_ASSERT(binParser.GetNextValue()==5);
	CPPUNIT_ASSERT(binParser.GetNextValue()==6);
	CPPUNIT_ASSERT(binParser.GetNextValue()==1);
	CPPUNIT_ASSERT(binParser.GetNextValue()==0);
	CPPUNIT_ASSERT(binParser.GetNextValue()==2);



}


}
}
