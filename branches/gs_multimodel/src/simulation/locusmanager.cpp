//
// C++ Implementation: locusmanager
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locusmanager.h"
#include "locus.h"
#include "locusxy.h"

namespace Simulation {
/******************************************************************************
 ***             Locus Array
 ******************************************************************************/



#ifdef CPPUNIT
CPPUNIT_TEST_SUITE_REGISTRATION(LocusManagerTest);

void WriteLocusFile(const char *filename) {
	ofstream file(filename);
	file<<"Test Locus Report\n6 Loci\nA bunch of labels\n";
	file<<"RL1-1\t0.5\t0.5\t0.0001\t100\tRL1-1\n";
	file<<"RL1-2\t0.5\t0.5\t0.0001\t200\tRL1-2\n";
	file<<"RL1-3\t0.3\t0.7\t0.0001\t300\tRL1-3\n";
	file<<"RL1-4\t0.3\t0.7\t0.0001\t400\tRL1-4\n";
	file<<"RL1-5\t0.8\t0.2\t0.01\t500\tRL1-5\n";
	file<<"RL1-6\t0.2\t0.8\t0.01\t600\tRL1-6";
	file.close();
}

void WriteLocusFile(const char *filename, int locusCount) {
	ofstream file(filename);
	file<<"Test Locus Report\n"<<locusCount<<" Loci\nA bunch of labels\n";
	float position = 0.0;
	for (int i=0; i<locusCount; i++) {
		int m = i%5+2;
		float maf = (float)m / 10.0;
		float dist = maf/100.0;
		position+=(maf * 100000.0);
		file<<"RL1-"<<i+1<<"\t"<<maf<<"\t"<<1.0-maf<<"\t"<<dist<<"\t"<<position<<"\t"<<"RL1-"<<i<<"\n";
	}
}

void WriteLocusXY(const char *filename) {
	ofstream file(filename);
	file<<"Test XY Locus\n15 Loci\nA bunch of labels\n";
	file<<"RL1-1\t0.5\t0.5\t0.00008\t0.01\t100\t100\tPAR\n";
	file<<"RL1-2\t0.2\t0.3\t0.00008\t0.01\t200\t200\tPAR\n";
	file<<"RL1-3\t0.3\t0.8\t0.00008\t0.10\t300\t350\tPAR\n";
	file<<"RL1-4\t0.2\t0.0\t0.00080\t0.10\t1100\t1100\tX_ONLY\n";
	file<<"RL1-5\t0.2\t0.0\t0.00008\t0.01\t1400\t1400\tX_ONLY\n";
	file<<"RL1-6\t0.0\t0.3\t0.00008\t0.01\t1500\t1500\tY_ONLY\n";
	file<<"RL1-7\t0.0\t0.3\t0.00001\t0.01\t1600\t1600\tY_HOMOLOG\n";
	file<<"RL1-8\t1.0\t0.3\t0.00001\t0.10\t2500\t2500\tY_HOMOLOG\n";
	file<<"RL1-9\t0.1\t0.0\t0.0008\t0.1\t2600\t2700\tX_HOMOLOG\n";
	file<<"RL1-10\t0.1\t0.0\t0.0008\t0.0001\t2700\t2800\tX_HOMOLOG\n";
	file<<"RL1-11\t0.2\t0.2\t0.00006\t0.0001\t3000\t3000\tPAR\n";
	file<<"RL1-12\t0.3\t0.3\t0.00005\t0.00008\t3080\t3080\tPAR\n";
	file<<"RL1-13\t0.4\t0.4\t0.00005\t0.00009\t3090\t3090\tPAR\n";
	file<<"RL1-14\t0.3\t0.3\t0.00003\t0.0001\t3130\t3130\tPAR\n";
	file<<"RL1-15\t0.25\t0.25\t0.00003\t0.0001\t3180\t3180\tPAR";
}
	
LocusManagerTest::LocusManagerTest() { }

LocusManagerTest::~LocusManagerTest() { }

void LocusManagerTest::setUp( ) { }

void LocusManagerTest::tearDown() { }


void LocusManagerTest::TestLocusXY() {
	string filename("test.loc");
	WriteLocusXY(filename.c_str());
	LocusManagerFileBased<LocusXY> fb("xytest", 1);
	fb.Load(filename.c_str());
	
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::LocusXY Count", 15, (int)fb.LocusCount());
	LocusXY *loc = fb[0];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Label(0)", strcmp(loc->GetLabel().c_str(), "RL1-1")==0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 1", 0.5, loc->Freq1X(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 2", 0.5, loc->Freq2X(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 1", 0.5, loc->Freq1Y(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 2", 0.5, loc->Freq2Y(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Map Units", 0.00008, loc->MapDistance(), 0.000001);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::LocusXY Location", 100, (int)loc->GetLocation());
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Label", strcmp("RL1-1", loc->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Description", strcmp("PAR", loc->GetDescription().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::PAR);

	loc = fb[4];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::X_Only);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 1", 0.2, loc->Freq1X(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 2", 0.8, loc->Freq2X(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 1", 0.0, loc->Freq1Y(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::LocusXY Freq 1", 1.0, loc->Freq2Y(), 0.001);
	
	loc=fb[5];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::Y_Only);
	loc=fb[6];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::Y_Homolog);
	loc=fb[8];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::X_Homolog);
	loc=fb[10];
	CPPUNIT_ASSERT_MESSAGE("Load::LocusXY Type", loc->type==LocusXY::PAR);
	
}

void LocusManagerTest::TestRegularLocus() {
	//Create a 5 locus array
	string filename("test.loc");
	WriteLocusFile(filename.c_str());
	LocusManagerFileBased<Locus> fb("test-project", 1);
	fb.Load(filename.c_str());	
	

	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::Locus Count", 6, (int)fb.LocusCount());
	Locus *loc = fb[0];
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Freq 1", 0.5, loc->Freq1(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Freq 2", 0.5, loc->Freq2(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Map Units", 0.0001, loc->MapDistance(), 0.000001);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::Locus Location", 100, (int)loc->GetLocation());
	CPPUNIT_ASSERT_MESSAGE("Load::Locus Label", strcmp("RL1-1", loc->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Load::Locus Description", strcmp("RL1-1", loc->GetDescription().c_str())==0);

	loc = fb[1];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::Locus Position", 200, (int)loc->GetLocation());
	
	loc = fb["RL1-3"];
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::Locus Position", 300, (int)loc->GetLocation());
	
	loc = fb.At(5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Freq 1", 0.2, loc->Freq1(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Freq 2", 0.8, loc->Freq2(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus Map Units", 0.01, loc->MapDistance(), 0.000001);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Load::Locus Location", 600, (int)loc->GetLocation());
	CPPUNIT_ASSERT_MESSAGE("Load::Locus Label", strcmp("RL1-6", loc->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Load::Locus Description", strcmp("RL1-6", loc->GetDescription().c_str())==0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Load::Locus MAF", 0.2, loc->GetMinAlleleFreq(), 0.000001);

	//Let's change a few allele frequencies
	loc->AssignFreq(0.1, 0.9);
	loc = fb.At(4);
	loc->SetFreq1(0.3);
	loc->SetFreq2(0.3);
	CPPUNIT_ASSERT_MESSAGE("!Locus Valid()", !loc->Valid());
	loc->SetFreq2(0.7);
	CPPUNIT_ASSERT_MESSAGE("Locus Valid()", loc->Valid());
	fb.Save(5, "testing");

	LocusManagerFileBased<Locus> newFb("test-project", 2);
	newFb.Load(filename.c_str());
	newFb.Refresh(5);

	loc = newFb.At(5);

	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Refresh::Locus Freq 1", 0.1, loc->Freq1(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Refresh::Locus Freq 2", 0.9, loc->Freq2(), 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Refresh::Locus Map Units", 0.01, loc->MapDistance(), 0.000001);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Refresh::Locus Location", 600, (int)loc->GetLocation());
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Refresh::Locus Genetic Location", 0.0104, loc->MapPosition(), 0.000001);
	CPPUNIT_ASSERT_MESSAGE("Refresh::Locus Label", strcmp("RL1-6", loc->GetLabel().c_str())==0);
	CPPUNIT_ASSERT_MESSAGE("Refresh::Locus Description", strcmp("RL1-6", loc->GetDescription().c_str())==0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Refresh::Locus MAF", 0.1, loc->GetMinAlleleFreq(), 0.000001);


	
}

#endif //CPPUNIT


}
