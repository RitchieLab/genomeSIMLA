//
// C++ Interface: peddatatest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_TESTPEDDATATEST_H
#define GENETICS_TESTPEDDATATEST_H

#include <cppunit/extensions/HelperMacros.h>

#include "genetics/familyrepository.h"
#include "utility/utility.h"

namespace Genetics {

namespace Test {

using namespace Utility;
using namespace Genetics;


/**
@brief Tests for basic functionality for the ped data structures

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PedDataTest: public CPPUNIT_NS::TestFixture {
	CPPUNIT_TEST_SUITE( PedDataTest );
	CPPUNIT_TEST( TestBasicInput );
	CPPUNIT_TEST( TestGenotypeData);
	CPPUNIT_TEST( TestErrorChecking );
	CPPUNIT_TEST( TestGenotypeData );
	CPPUNIT_TEST( TestMedelianError );
	CPPUNIT_TEST( TestMendelianErrorClearingSibs );
	CPPUNIT_TEST( TestMendelianErrorRemoves );
	CPPUNIT_TEST_SUITE_END();

public:
    PedDataTest();

    ~PedDataTest();
	
	void setUp();
	
	void TestGenotypeData();
	void TestBasicInput();
	void TestErrorChecking();
	void TestFamilyArrangements();
	void TestMedelianError();
	void TestTrioFix();
	void TestMendelianErrorRemoves();
	void TestMendelianErrorClearingSibs();
	void tearDown();

protected:
	FamilyRepository *familyRepo;

};



}

}

#endif
