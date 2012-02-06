//
// C++ Implementation: allelesource
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "allelesource.h"
#include "locus.h"
#include "locusxy.h"
#include "locusmanager.h"
#include <string>

using namespace std;
namespace Simulation {

bool DoSuspend = false;


#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(AlleleSourceTest);


AlleleSourceTest::AlleleSourceTest() { }

AlleleSourceTest::~AlleleSourceTest() { }

void AlleleSourceTest::setUp() { }

void AlleleSourceTest::tearDown() { }

void AlleleSourceTest::TestCrossOverXY() {
	string filename("20locXY.loc");
	WriteLocusXY(filename.c_str());
	LocusManagerFileBased<LocusXY> fb("XY", 1);
	fb.Load(filename.c_str());

	AlleleSourceSingle<LocusXY> *allZero=new AlleleSourceSingle<LocusXY>(&fb, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle<XY>::Locus Count", 15, (int)allZero->GetLocusCount());
	AlleleSourceSingle<LocusXY> *allOnes=new AlleleSourceSingle<LocusXY>(&fb, 2);
	for (int i=0; i<15; i++) {
		allZero->At(i, 0);
		allOnes->At(i, 1);
	}

	Utility::RBTree<int, AlleleSource<LocusXY>*> pool;
	pool.Add(allZero->GetSourceID(), allZero);
	pool.Add(allOnes->GetSourceID(), allOnes);
	
	Utility::RBTree<float, AlleleSource<LocusXY>*> sources;
	allZero->GetSource(sources, 0, 14);
	
	CPPUNIT_ASSERT_EQUAL_MESSAGE("GetSource()", 1, sources.GetCount());
	
	vector<size_t> xoPoints;
	xoPoints.push_back(2);
	xoPoints.push_back(3);
	xoPoints.push_back(14);
	AlleleSource<LocusXY> *newCopy = allZero->Cross(allOnes, xoPoints);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(3));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(4));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(5));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(13));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(14));

	delete newCopy;
	delete allZero;
	delete allOnes;
	
}
void AlleleSourceTest::TestCrossOver() {
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
	for (int i=0; i<20; i++) {
		allZero->At(i, 0);
		allOnes->At(i, 1);
	}

	//Create a pool of our two chromosomes
	Utility::RBTree<int, AlleleSource<Locus>*> pool;
	pool.Add(allZero->GetSourceID(), allZero);
	pool.Add(allOnes->GetSourceID(), allOnes);

	//Check out the GetSource Stuff
	Utility::RBTree<float, AlleleSource<Locus>* > sources;
	allZero->GetSource(sources, 0, 19);
	//There should only be one item
	CPPUNIT_ASSERT_EQUAL_MESSAGE("GetSource()", 1, sources.GetCount());
	
	sources.Clear();
	
	vector<size_t> xoPoints;
	xoPoints.push_back(5);
	xoPoints.push_back(10);
	xoPoints.push_back(15);
	AlleleSource<Locus> *newCopy = allZero->Cross(allOnes, xoPoints);
	//Let's check on the loci...it seems to be getting mixed up here for some reason
	LocusManager<Locus> * loci = allZero->GetLoci();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Locus after cross over", 20, (int)allZero->GetLocusCount());
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-1", loci->At(0)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-2", loci->At(1)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-3", loci->At(2)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-4", loci->At(3)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-5", loci->At(4)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-6", loci->At(5)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-7", loci->At(6)->GetLabel().c_str())==0);
	loci = newCopy->GetLoci();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Locus after cross over", 20, (int)newCopy->GetLocusCount());
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-1", loci->At(0)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-2", loci->At(1)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-3", loci->At(2)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-4", loci->At(3)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-5", loci->At(4)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-6", loci->At(5)->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Locus after cross over", strcmp("RL1-7", loci->At(6)->GetLabel().c_str())==0);

	newCopy->GetSource(sources, 0, 19);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("GetSource()", 4, sources.GetCount());		
	CPPUNIT_ASSERT_EQUAL_MESSAGE("GetSource()->First", allZero->GetSourceID(), sources.GetFirst()->GetData()->GetSourceID());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("GetSource()->Last", allOnes->GetSourceID(), sources.GetLast()->GetData()->GetSourceID());
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(3));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(4));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(5));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, newCopy->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, newCopy->At(19));

	

	xoPoints.clear();
	xoPoints.push_back(0);
	xoPoints.push_back(2);
	AlleleSource<Locus> *sibling = allZero->Cross(allOnes, xoPoints);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, sibling->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, sibling->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, sibling->At(2));
	


	filename = "tempfile.bin";
	ofstream file(filename.c_str(), ios::trunc|ios::binary);
	newCopy->WriteBinary(file);
	sibling->WriteBinary(file);
	newCopy->WriteCompleteBinary(file, -1);
	sibling->WriteCompleteBinary(file, -1);
	file.close();
	
	ifstream infile(filename.c_str(), ios::in|ios::binary);
	AlleleSource<Locus> *twin1 = new AlleleSourceSingle<Locus>(&fb, 3);
	AlleleSource<Locus> *twin2 = new AlleleSourceSingle<Locus>(&fb, 4);
	AlleleSource<Locus> *twin3 = new AlleleSourceInh<Locus>(&fb, 5);
	AlleleSource<Locus> *twin4 = new AlleleSourceInh<Locus>(&fb, 6);
	twin3->ReadBinary(infile, &pool);
	twin4->ReadBinary(infile, &pool);
	twin1->ReadBinary(infile, &pool);
	twin2->ReadBinary(infile, &pool);	
	infile.close();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin1->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin1->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin1->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin1->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin1->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin1->At(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin1->At(19));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin2->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin2->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin2->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin3->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin3->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin3->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin3->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin3->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin3->At(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin3->At(19));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin4->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, twin4->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, twin4->At(2));

	
	
	file.open(filename.c_str(), ios::in|ios::binary);

	AlleleSource<Locus> *copyOfACopy = newCopy->Clone();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, copyOfACopy->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, copyOfACopy->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, copyOfACopy->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, copyOfACopy->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, copyOfACopy->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, copyOfACopy->At(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, copyOfACopy->At(19));
	
	AlleleSource<Locus> *realized = newCopy->Realize(7);

	delete allZero;
	delete allOnes;
	//OK, we've deleted the source structures...this should only be possible, if realize worked
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, realized->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, realized->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, realized->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, realized->At(9));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 0, realized->At(10));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, realized->At(15));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::XO", 1, realized->At(19));
	

	delete newCopy;
	delete copyOfACopy;
	delete sibling;
	delete twin1;
	delete twin2;
	delete twin3;
	delete twin4;
	delete realized;
}
void AlleleSourceTest::TestInitialization() {
	string filename("test.loc");
	WriteLocusFile(filename.c_str());
	LocusManagerFileBased<Locus> fb("test", 1);
	fb.Load(filename.c_str());

	AlleleSourceSingle<Locus> *allSrc = new AlleleSourceSingle<Locus>(&fb, 1);
	allSrc->At(0, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, allSrc->At(0));
	allSrc->At(1, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, allSrc->At(1));
	allSrc->At(2, 1);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, allSrc->At(2));
	allSrc->At(3, 0);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, allSrc->At(3));
	allSrc->At(4, 0);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, allSrc->At(4));
	allSrc->At(5, 0);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, allSrc->At(5));

	AlleleSource<Locus> *copy = allSrc->Clone();
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, copy->At(0));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, copy->At(1));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 1, copy->At(2));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, copy->At(3));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, copy->At(4));
	CPPUNIT_ASSERT_EQUAL_MESSAGE("AlleleSourceSingle::at", 0, copy->At(5));
	
	delete copy;
	delete allSrc;
	
}

#endif //CPPUNIT

}
