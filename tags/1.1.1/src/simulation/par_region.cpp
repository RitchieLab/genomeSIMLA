//
// C++ Implementation: par_region
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "par_region.h"
#include "locusxy.h"
namespace Simulation {
#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(PARTest);

PARTest::PARTest() { }

PARTest::~PARTest() { }

void PARTest::setUp() { }

void PARTest::tearDown() { }


void PARTest::TestPAR() {
	PAR_Region<LocusXY> par;
	LocusXY l1(0.0, 1, 1, 0);
	l1.SetMapPositionY(0.0);
	LocusXY l2(0.00001, 1, 2, 10);
	l2.SetMapPositionY(0.00001);
	LocusXY l3(0.00001, 1, 3, 20);
	l3.SetMapPositionY(0.00002);
	LocusXY l4(0.000001, 1, 4, 21);
	l4.SetMapPositionY(0.000021);
	LocusXY l5(0.000001, 1, 10, 1540);
	l5.SetMapPositionY(0.100000);
	LocusXY l6(0.000005, 1, 11, 1545);
	l6.SetMapPositionY(0.100005);

	par.InsertRegion(l1, l2);
	par.InsertRegion(l2, l3);
	par.InsertRegion(l3, l4);
	par.InsertRegion(l5, l6);
	
}

#endif //CPPUNIT
}
