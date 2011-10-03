//
// C++ Implementation: cpool
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "cpool.h"
#include "cpoolxy.h"
#include "locus.h"
#include "locusmanager.h"
#include "utility/random.h"
#include "lineargrowth.h"


namespace Simulation {

unsigned int targetPop					=0;
unsigned int numberOfBlocksToReport		=15;	///<Control the number of blocks written to the final report
unsigned int highBufferSize 			=50;	///<The zoomed out buffer around the target block
unsigned int detailsBufferSize			=5;		///<The smaller buffer 

#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(CPoolTest);


CPoolTest::CPoolTest() { }

CPoolTest::~CPoolTest() { }

void CPoolTest::setUp() { }

void CPoolTest::tearDown() { }


void CPoolTest::TestInitialization() {
	Utility::Random rnd;
	string filename("20loc.loc");
	WriteLocusFile(filename.c_str(), 20);
	LocusManagerFileBased<Locus> fb("test", 1);
	fb.Load(filename.c_str());

	//I have no idea how to test randomized stuff. So, I'm going to just exercise 
	//the functions and look at things I can test for. Verifying things like
	//the genotypes inside a given individual, though, is going to be hard to do
	CPool<Locus> pool(1, "test_pool", "Test Pool");
	pool.Init(&fb, rnd);
	pool.BuildInitialPopulation(rnd, 100);
	pool.Save();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", 100, (int)pool.GetPopulationSize());
	//OK, let's try opening the pool as another and make sure that they have the same
	//stuff inside
	CPool<Locus> pool2(2, "test_pool", "Test Pool");
	pool2.Open(&fb, rnd);
	
	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", 100, (int)pool2.GetPopulationSize());
	AlleleSource<Locus> *chrom1 = pool.At(0);
	AlleleSource<Locus> *chrom2 = pool2.At(0);

	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", 20, (int)chrom2->GetLocusCount());
	int countZeros=0;
	for (int i=0; i<20; i++) {
		CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", chrom1->At(i), chrom2->At(i));
		countZeros+=(chrom1->At(i)==1);
	}
	CPPUNIT_ASSERT_MESSAGE("num zeros", (countZeros > 0 && countZeros < 20));

	//These should be different once we advance into the future
	chrom1 = pool.At(0);
		int a = chrom1->At(0);
	PopulationGrowth::GrowthRate::minPoolSize = 100;
	//Grow the population
	PopulationGrowth::LinearGrowth growth(100, 5, 0.0);				///<Basic growth, 5/generation
	pool.AdvanceGenerations(10, &growth, rnd, 1);
	

	CPPUNIT_ASSERT_EQUAL_MESSAGE("cpool::open", 150, (int)pool.GetPopulationSize());
	chrom2 = pool.At(0);
	
	countZeros = 0;
	int overlap =0;
	for (int i=0; i<20; i++) {
		int a = chrom1->At(i);
		int b = chrom2->At(i);
		overlap += chrom1->At(i) == chrom2->At(i);
		countZeros += chrom2->At(i) == 1;
	}
	CPPUNIT_ASSERT_MESSAGE("cpool::advance generation", (overlap > 0 && overlap < 20));
	CPPUNIT_ASSERT_MESSAGE("cpool::advance generation", (countZeros > 0 && countZeros < 20));
		
	pool.Suspend();
	pool.Wake();

	chrom2 = pool.At(0);

	int zeros=0, overlapping=0;
	for (int i=0; i<20; i++) {
		overlapping += chrom1->At(i) == chrom2->At(i);
		zeros += chrom2->At(i) == 1;
	}

}


#endif

}
