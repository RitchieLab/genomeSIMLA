//
// C++ Interface: testpdttstat
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESE_TESTTESTPDTTSTAT_H
#define ESE_TESTTESTPDTTSTAT_H

#include "genetics/gtfileparserbuffered.h"
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "gtlineparsermdrpdt.h"
#include "evalbalancedaccuracypdt.h"
#include "genetics/snprepository.h"

namespace ESE {

namespace Test {

using namespace std;
using namespace Genetics::Parser;

/**
@brief Load several different data sets and verify the T-Statistic matches what we expect

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class TestPdtTStat : public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( TestPdtTStat );
	CPPUNIT_TEST( TestModels2PCx31 );
	CPPUNIT_TEST( TestModels1x40 );
	CPPUNIT_TEST( TestModels1x35 );
	CPPUNIT_TEST( TestModels2PCx12 );
	CPPUNIT_TEST( TestModelsEpiAAU1x76 );
	CPPUNIT_TEST( TestModelsEpiAAU1x100 );
	CPPUNIT_TEST_SUITE_END();
public:
    TestPdtTStat();

    ~TestPdtTStat();
	void TestModels1x35();
	void TestModels1x40();
	void TestModels2PCx31();
	void TestModels2PCx12();
	void TestModelsEpiAAU1x76();
	void TestModelsEpiAAU1x100();

protected:
	void SetupRepository(const char *filename, bool doCreateVSibs = false);
	void CleanupRepository();
	void VerifyTScores(string models[], float tValues[], uint modelCount );
	void VerifyGenotypeCounts(const char *modelID, uint affectedCount, uint unaffectedCount, uint expValues[] );
	FamilyRepository *families;
	GtLineParserMdrPdt *pdtParser;
	GtFileParserBuffered *fileParser;
	CaseControlStatus status;
	SnpRepository *snps;
	EvalBalancedAccuracyPDT *evalPDT;
	FoldType pgStatus;


};

}

}

#endif
