//
// C++ Implementation: cpair.cpp
//
// Description: 
//
//
// Author:  <Eric Torstenson>, Marylyn Ritchie (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "cpair.h"

namespace Simulation {

#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(CPairTest);


CPairTest::CPairTest() { }

CPairTest::~CPairTest() { }

void CPairTest::setUp() { }

void CPairTest::tearDown() { }


void CPairTest::Test() {
	string filename("20loc.loc");
	WriteLocusFile(filename.c_str(), 20);
	LocusManagerFileBased<Locus> fb("test", 1);
	fb.Load(filename.c_str());

	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-1", fb.At(0)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-2", fb.At(1)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-3", fb.At(2)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-4", fb.At(3)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-5", fb.At(4)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-6", fb.At(5)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-7", fb.At(6)->GetLabel().c_str())==0);
	
	AlleleSourceSingle<Locus> *allZero=new AlleleSourceSingle<Locus>(&fb, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::Locus Count", 20, (int)allZero->GetLocusCount());
	AlleleSourceSingle<Locus> *allOnes=new AlleleSourceSingle<Locus>(&fb, 2);
	AlleleSourceSingle<Locus> *halfNHalf=new AlleleSourceSingle<Locus>(&fb, 3);
	for (int i=0; i<20; i++) {
		allZero->At(i, 0);
		allOnes->At(i, 1);
		halfNHalf->At(i, i<10);
	}
	
	vector<size_t> xoPoints;
	xoPoints.push_back(5);
	xoPoints.push_back(10);
	xoPoints.push_back(15);
	AlleleSource<Locus> *newCopy = allZero->Cross(allOnes, xoPoints);
	xoPoints.clear();
	xoPoints.push_back(10);
	AlleleSource<Locus> *newHalf = allZero->Cross(allOnes, xoPoints);
	
	//Normally, you don't want to clone these for pairs, but I'm reusing them, so I have to
	CPair<Locus> pair1(allZero->Clone(), halfNHalf->Clone(), &fb);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair1.GetGenotype(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair1.GetGenotype(19));
	
	CPair<Locus> pair2(allOnes->Clone(), halfNHalf->Clone(), &fb);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair2.GetGenotype(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair2.GetGenotype(19));

	CPair<Locus> pair3(newHalf->Clone(), newCopy->Clone(), &fb);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(4));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(5));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(14));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair3.GetGenotype(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair3.GetGenotype(19));

	pair3.SetMissing(4, true);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype (missing)", -1, pair3.GetGenotype(4));

	delete newHalf;
	delete newCopy;
	delete allZero;
	delete allOnes;
	delete halfNHalf;
	
}

void CPairTest::TestXY() {
	string filename("20locXY.loc");
	WriteLocusXY(filename.c_str());
	LocusManagerFileBased<LocusXY> fb("XY", 1);
	fb.Load(filename.c_str());
	
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-1", fb.At(0)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-2", fb.At(1)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-3", fb.At(2)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-4", fb.At(3)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-5", fb.At(4)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-6", fb.At(5)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-7", fb.At(6)->GetLabel().c_str())==0);


	AlleleSourceSingle<LocusXY> *allZero=new AlleleSourceSingle<LocusXY>(&fb, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle<XY>::Locus Count", 15, (int)allZero->GetLocusCount());
	AlleleSourceSingle<LocusXY> *allOnes=new AlleleSourceSingle<LocusXY>(&fb, 2);
	AlleleSourceSingle<LocusXY> *halfNHalf=new AlleleSourceSingle<LocusXY>(&fb, 3);
	for (int i=0; i<15; i++) {
		allZero->At(i, 0);
		allOnes->At(i, 1);
		halfNHalf->At(i, (i<8));
	}


	
	vector<size_t> xoPoints;
	xoPoints.push_back(4);
	xoPoints.push_back(8);
	xoPoints.push_back(12);

	AlleleSource<LocusXY> *newCopy = allZero->Cross(allOnes, xoPoints);

	xoPoints.clear();
	xoPoints.push_back(7);
	AlleleSource<LocusXY> *newHalf = allZero->Cross(allOnes, xoPoints);
	
	//Normally, you don't want to clone these for pairs, but I'm reusing them, so I have to
	CPairXY pair1(allZero->Clone(), halfNHalf->Clone(), &fb, true);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::IsXX", true, pair1.IsXX());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair1.GetGenotype(7));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair1.GetGenotype(8));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair1.GetGenotype(14));
	
	CPairXY pair2(allOnes->Clone(), halfNHalf->Clone(), &fb, false);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::IsXX", false, pair2.IsXX());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair2.GetGenotype(7));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair2.GetGenotype(8));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair2.GetGenotype(14));

	CPairXY pair3(newHalf->Clone(), newCopy->Clone(), &fb, true);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::IsXX", true, pair3.IsXX());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 0, pair3.GetGenotype(3));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(4));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair3.GetGenotype(7));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(8));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 1, pair3.GetGenotype(11));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair3.GetGenotype(12));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("CPairTest::GetGenotype", 2, pair3.GetGenotype(14));

	delete newHalf;
	delete newCopy;
	delete allZero;
	delete allOnes;
	delete halfNHalf;
	
}
#endif


}

