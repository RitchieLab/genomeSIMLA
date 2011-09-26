//
// C++ Implementation: testcasecontrolstatus
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "testcasecontrolstatus.h"

namespace Utility {

namespace Test {



using namespace std;
CPPUNIT_TEST_SUITE_REGISTRATION( TestCaseControlStatus );

void TestCaseControlStatus::TestCopyConstructor() {
	BitSetType a(10, false);
	BitSetType u(10, false);
	BitSetType t(10, true);
	for (uint i=0; i<5; i++) {
		a[i]=true;
		t[i]=true;
	}

	for (uint i=5; i<10; i++) {
		u[i]=true;
		t[i]=true;
	}

	CaseControlStatus aa(a, u);
	CaseControlStatus cc(aa);

	CPPUNIT_ASSERT(cc.affected.count() == 5);
	CPPUNIT_ASSERT(cc.unaffected.count() == 5);
	cout<<"\nTestCaseControlStatus::t: "<<t<<"\n";
	cout<<"TestCaseControlStatus::cc.total: "<<cc.total<<"\n";
	CPPUNIT_ASSERT(cc.total == t);
	
	CaseControlStatus dd=cc.MakeRandomCopy();
	cout<<"\nTestCaseControlStatus::dd.affected:   "<<dd.affected<<"\n";
	cout<<  "TestCaseControlStatus::dd.unaffected: "<<dd.unaffected<<"\n";
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 5);
	CPPUNIT_ASSERT(dd.affected != a);
	CPPUNIT_ASSERT(dd.unaffected != u);
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 5);
	CPPUNIT_ASSERT((dd.affected & dd.unaffected).count()==0);	
}

void TestCaseControlStatus::TestBalanced() {
	BitSetType a(10, false);
	BitSetType u(10, false);
	BitSetType t(10, true);
	for (uint i=0; i<5; i++) {
		a[i]=true;
		t[i]=true;
	}

	for (uint i=5; i<10; i++) {
		u[i]=true;
		t[i]=true;
	}

	CaseControlStatus cc(a, u);
	CPPUNIT_ASSERT(cc.affected.count() == 5);
	CPPUNIT_ASSERT(cc.unaffected.count() == 5);
	cout<<"\nTestCaseControlStatus::t: "<<t<<"\n";
	cout<<"TestCaseControlStatus::cc.total: "<<cc.total<<"\n";
	CPPUNIT_ASSERT(cc.total == t);
	
	CaseControlStatus dd=cc.MakeRandomCopy();
	cout<<"\nTestCaseControlStatus::dd.affected:   "<<dd.affected<<"\n";
	cout<<  "TestCaseControlStatus::dd.unaffected: "<<dd.unaffected<<"\n";
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 5);
	CPPUNIT_ASSERT(dd.affected != a);
	CPPUNIT_ASSERT(dd.unaffected != u);
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 5);
	CPPUNIT_ASSERT((dd.affected & dd.unaffected).count()==0);
	 
}

void TestCaseControlStatus::TestImperfect() {
	BitSetType a(15, false);
	BitSetType u(15, false);
	BitSetType t(15, false);
	t[5]=false;
	t[14]=false;
	//5 affected
	for (uint i=0; i<5; i++) { 
		a[i]=true;
		t[i]=true;
	}

	//8 unaffected
	for (uint i=6; i<14; i++)	{
		u[i]=true;
		t[i]=true;
	}

	CaseControlStatus cc(a, u);
	CPPUNIT_ASSERT(cc.affected.count() == 5);
	CPPUNIT_ASSERT(cc.unaffected.count() == 8);
	CPPUNIT_ASSERT(cc.total == t);
	
	CaseControlStatus dd=cc.MakeRandomCopy();
	cout<<"\nTestCaseControlStatus::dd.affected:   "<<dd.affected<<"\n";
	cout<<  "TestCaseControlStatus::dd.unaffected: "<<dd.unaffected<<"\n";
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 8);
	CPPUNIT_ASSERT(dd.affected != a);
	CPPUNIT_ASSERT(dd.unaffected != u);
	CPPUNIT_ASSERT(dd.affected.count() == 5);
	CPPUNIT_ASSERT(dd.unaffected.count() == 8);

	//Now, we have a couple we know CAN'T be true
	CPPUNIT_ASSERT(!dd.affected[5]);
	CPPUNIT_ASSERT(!dd.affected[15]);
	CPPUNIT_ASSERT(!dd.unaffected[5]);
	CPPUNIT_ASSERT(!dd.unaffected[15]);	
	CPPUNIT_ASSERT((dd.affected & dd.unaffected).count()==0);

}

void TestCaseControlStatus::TestSetStatus() {
	CaseControlStatus a(0);
	a.SetStatus(15, true);
	a.SetStatus(13, false);
	a.SetStatus(1, true);
	a.SetStatus(2, false);
	a.SetStatus(3, true);
	a.SetStatus(5, true);
	a.SetStatus(10, false);
	
	BitSetType aff(16, false);
	BitSetType unaff(16, false);
	unaff[13]=true;
	unaff[2]=true;
	unaff[10] = true;
	aff[1]=true;
	aff[3]=true;
	aff[5]=true;
	aff[15]=true;
	

	CPPUNIT_ASSERT(a.affected.size() == 16);
	CPPUNIT_ASSERT(a.affected.count() == 4);
	CPPUNIT_ASSERT(a.unaffected.size() == 16);
	CPPUNIT_ASSERT(a.unaffected.count() == 3);
	CPPUNIT_ASSERT(a.total.count() == 7);
	CPPUNIT_ASSERT(a.total.size() == 16);
	CPPUNIT_ASSERT(a.total == (aff | unaff));
	CPPUNIT_ASSERT(a.unaffected == unaff);
	CPPUNIT_ASSERT(a.affected == aff);
	
	cout<<"TEsting.";
}
void TestCaseControlStatus::TestAppendStatus() {
	BitSetType a(10, false);
	BitSetType c(20, false);
	BitSetType u(10, false);
	BitSetType cu(20, false);
	for (uint i=0; i<5; i++) {
		a[i]=true;
		c[i+10]=true;
	}

	for (uint i=5; i<10; i++) {
		u[i]=true;
		cu[i+10]=true;
	}

	CaseControlStatus aa(a, u);
	CaseControlStatus bb(c, cu);
	aa.AppendStatus(bb);

	a.resize(20);
	u.resize(20);
	for (uint i=10; i<15; i++) {
		a[i]=true;
		u[i+5]=true;
	}	

	CPPUNIT_ASSERT(aa.affected.size() == 20);
	CPPUNIT_ASSERT(aa.affected.count() == 10);
	CPPUNIT_ASSERT(aa.unaffected.count() == 10);
	CPPUNIT_ASSERT(aa.affected == a);
	cout<<"\nAA : "<<aa.affected<<"\n";
	cout<<"\na  : "<<a<<"\n";
	cout<<"\nAA : "<<aa.unaffected<<"\n";
	cout<<"u  : "<<u<<"\n";
	CPPUNIT_ASSERT(aa.unaffected == u);
	cout<<"Testing 1 2 3\n";

}

}

}
