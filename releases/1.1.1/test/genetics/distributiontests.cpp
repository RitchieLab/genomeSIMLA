//
// C++ Implementation: distributiontests
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "distributiontests.h"
#include "genetics/ptestdistribution.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Genetics::Test::DistributionTests);


namespace Genetics {

namespace Test {


DistributionTests::DistributionTests() {
}


DistributionTests::~DistributionTests()
{
}
/*
*       0      3.742
*       1      3.642
*       2      3.450
*       3      3.329
*       4      3.265
*       5      3.220
*       6      3.014
*       7      3.007
*       8      2.952
*       9      2.917


*       0      4.604
*       1      4.492
*       2      4.483

*       3      4.471
*       4      4.457
*       5      4.391
*       6      4.391
*       7      4.337
*       8      4.298
*       9      4.192
*/


void DistributionTests::TestNTestDistribution() {
	Distribution::NTest dist(2, 13, "N-Test Distribution Test");
	float fakeScores[2][10] = { {3.014,2.917,3.329,3.265,3.642,3.220,3.742,3.450,3.007,2.952},
							   {4.391,4.337,4.471,4.457,4.604,4.483,4.298,4.492,4.391,4.192} };
	char id[6];
	//Let's add a few numbers to the distribution
	for (uint inner=0; inner<10; inner++) {
		sprintf(id, "%d", inner);
		dist.AppendTest( fakeScores[0][inner], 0, id, inner);
		sprintf(id, "2x%d", inner);
		dist.AppendTest( fakeScores[1][inner], 1, id, inner);
	}
	

	//We can't get pvalues from unsorted, so we might as well sort it
	dist.Sort();
	cout<<"\ndist.GetPValue(1, 3.224)="<<dist.GetPValue( 3.224, 1)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 3.224, 1), 0.6000, 0.001);

	cout<<"dist.GetPValue(1, 3.8)="<<dist.GetPValue( 3.8, 1)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 3.8, 1), 0.1000, 0.001);

	cout<<"dist.GetPValue(1, 2.1)="<<dist.GetPValue( 2.1, 1)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 2.1, 1), 1.000, 0.001);

	cout<<"dist.GetPValue(1, 2.1)="<<dist.GetPValue( 2.1, 1)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 2.952, 1), 0.9000, 0.001);

	cout<<"dist.GetPValue(2, 4.399)="<<dist.GetPValue( 4.399, 2)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 4.399, 2), 0.6000, 0.001);

	cout<<"dist.GetPValue(2, 4.65)="<<dist.GetPValue( 4.65, 2)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 4.650, 2), 0.1000, 0.001);

	cout<<"dist.GetPValue(2, 4.480)="<<dist.GetPValue( 4.480, 2)<<"\n";
	CPPUNIT_ASSERT_DOUBLES_EQUAL(dist.GetPValue( 4.480, 2), 0.4000, 0.001);


}

}

}
