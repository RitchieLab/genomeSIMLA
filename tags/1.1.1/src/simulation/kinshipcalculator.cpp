//
// C++ Implementation: kinshipcalculator
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kinshipcalculator.h"

namespace Simulation {


#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(KinshipCalculatorTest);

KinshipCalculatorTest::KinshipCalculatorTest() {	}


KinshipCalculatorTest::~KinshipCalculatorTest(){	}
void KinshipCalculatorTest::tearDown() {

}


void KinshipCalculatorTest::setUp() {
	left.Insert (0,  "P");
	left.Insert (5,  "M");
	left.Insert (8,  "P");
	left.Insert (12, "M");
	left.Insert (13, "P");
	left.Insert (15, "M");
	left.Insert (19, "P");
	left.Insert (21, "X");
	
	right.Insert(0,  "M");
	right.Insert(3,  "P");
	right.Insert(9,  "M");
	right.Insert(14, "P");
	right.Insert(21, "X");
}


void KinshipCalculatorTest::TestSegment() {
	KinshipCalculator<string>::Segment seg1(0,5);
	KinshipCalculator<string>::Segment seg2(3, 8);
	KinshipCalculator<string>::Segment seg3(1, 5);
	int direction;
	CPPUNIT_ASSERT_EQUAL(0, seg1.start);
	CPPUNIT_ASSERT_EQUAL(8, seg2.stop);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Overlap", 3, seg1.EvaluateOverlap(seg2, direction));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Direction", -1, direction);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Overlap", 6, seg1.EvaluateOverlap(seg1, direction));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Direction", -1, direction);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Overlap", 3, seg2.EvaluateOverlap(seg1, direction));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Direction", 1, direction);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Overlap", 5, seg1.EvaluateOverlap(seg3, direction));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Direction", -1, direction);

}

void KinshipCalculatorTest::TestSegmentation() {
	KinshipCalculator<string> k;
	KinshipCalculator<string>::SegmentMap lSeg;

	k.SegmentTree(left, lSeg);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("TestSegmentation (Left)", 2, (int)lSeg.size());
	KinshipCalculator<string>::SegmentVector &seg = lSeg["P"];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("TestSegmentation Contents (Left)", 4, (int)seg.size());
	CPPUNIT_ASSERT_MESSAGE("Segment[0]", seg[0]==KinshipCalculator<string>::Segment(0,4));
	CPPUNIT_ASSERT_MESSAGE("Segment[1]", seg[1]==KinshipCalculator<string>::Segment(8,11));
	CPPUNIT_ASSERT_MESSAGE("Segment[2]", seg[2]==KinshipCalculator<string>::Segment(13,14));
	CPPUNIT_ASSERT_MESSAGE("Segment[3]", seg[3]==KinshipCalculator<string>::Segment(19,20));

	seg =lSeg["M"];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("TestSegmentation Contents (Left)", 3, (int)seg.size());
	CPPUNIT_ASSERT_MESSAGE("Segment[0]", seg[0]==KinshipCalculator<string>::Segment(5,7));
	CPPUNIT_ASSERT_MESSAGE("Segment[1]", seg[1]==KinshipCalculator<string>::Segment(12,12));
	CPPUNIT_ASSERT_MESSAGE("Segment[2]", seg[2]==KinshipCalculator<string>::Segment(15,18));

	KinshipCalculator<string>::SegmentMap rSeg;
	k.SegmentTree(right, rSeg);

	seg = rSeg["P"];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Contents (Right)", 2, (int)seg.size());
	CPPUNIT_ASSERT_MESSAGE("Segment[0]", seg[0]==KinshipCalculator<string>::Segment(3,8));
	CPPUNIT_ASSERT_MESSAGE("Segment[1]", seg[1]==KinshipCalculator<string>::Segment(14,20));
	
	seg = rSeg["M"];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Contents (Right)", 2, (int)seg.size());
	CPPUNIT_ASSERT_MESSAGE("Segment[0]", seg[0]==KinshipCalculator<string>::Segment(0,2));
	CPPUNIT_ASSERT_MESSAGE("Segment[1]", seg[1]==KinshipCalculator<string>::Segment(9,13));
	
	
}
	
void KinshipCalculatorTest::TestKinship() {
	//

	

	KinshipCalculator<string> k;
	int kinship = k.EvaluateKinship(left, right);
	
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Kinship", 7, kinship);

	//I'm not sure if we want to allow 0 or not. Right now, we are....
	CPPUNIT_ASSERT_EQUAL_MESSAGE("kinship", 21, k.EvaluateKinship(left, left));
}

#endif

}
