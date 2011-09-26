//
// C++ Implementation: snppooltest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snppooltest.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Genetics::Test::SnpPoolTest );

namespace Genetics {

namespace Test {



SnpPoolTest::SnpPoolTest()
{
}


SnpPoolTest::~SnpPoolTest()
{
	pool->Release();
}
void SnpPoolTest::setUp() {
	SnpPool::Initialize(5, 5, 4);
	pool = SnpPool::Instance();
}

void SnpPoolTest::tearDown() {
	pool->Release();
}


void SnpPoolTest::TestAcquisitionAndRelease() {
	int snpCount=pool->GetUnusedCount();
	SnpAligned *snp = pool->GetSnp(1);
	
	CPPUNIT_ASSERT_MESSAGE("GetSnp() Should always return a valid pointer", snp != NULL);
	CPPUNIT_ASSERT_MESSAGE("There should be 1 less snp available", pool->GetUnusedCount() == snpCount-1);
	pool->ReleaseSnp(snp);
	CPPUNIT_ASSERT_MESSAGE("There should be the same number as when we started", pool->GetUnusedCount() == snpCount);

}
}

}
