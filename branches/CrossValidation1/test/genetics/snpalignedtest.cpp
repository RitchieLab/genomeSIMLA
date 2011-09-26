//
// C++ Implementation: snpalignedtest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snpalignedtest.h"
#include "genetics/snpaligned.h"
#include <iostream>
#include "utility/utility.h"


using namespace std;
using namespace Utility;
CPPUNIT_TEST_SUITE_REGISTRATION( Genetics::Test::SnpAlignedTest );

namespace Genetics {
namespace Test {

SnpAlignedTest::SnpAlignedTest()
{
}


SnpAlignedTest::~SnpAlignedTest()
{
	pool->Release();
}
void SnpAlignedTest::setUp() {
	SnpPool::Initialize(5, 5, 4);
	pool = SnpPool::Instance();
}

void SnpAlignedTest::tearDown() {
	pool->Release();
}

void SnpAlignedTest::TestDescriptors() {
	//The mask is used to determine the affected individuals. As we are setting upt the 
	//bitvectors, they read right->left.
	BitSetType mask(16, 255ul);			//0000000011111111
	char c1[]="2213122132222313";		//GA
	char c2[]="3222233211111111";		//TA
	char c3[]="2321312322231321";		//GC

	SnpAligned *s1=pool->GetSnp(4);
	SnpAligned *s2=pool->GetSnp(4);
	SnpAligned *s3=pool->GetSnp(4);
//	SnpAligned *s4=NULL;		//1x2
//	SnpAligned *s5=NULL;		//3x2
	/**
	 * <UL><LI>000 - M
	 * <LI>001 - A
	 * <LI>010 - C
 	 * <LI>011 - G
	 * <LI>100 - T
	 * </UL>
	 */

	s1->ImportSnp(0, 1, c1);		//0,0,0
	s1->SetBinaryGenotypes(25);
	s2->ImportSnp(1, 2, c2);		//-8,3,5
	s2->SetBinaryGenotypes(33);
	s3->ImportSnp(2, 3, c3);		//0,0,0
	s3->SetBinaryGenotypes(26);

	CPPUNIT_ASSERT(strcmp(s1->GetTxtDescriptor().c_str(), "1GA")==0);
	CPPUNIT_ASSERT(strcmp(s2->GetTxtDescriptor().c_str(), "2TA")==0);
	CPPUNIT_ASSERT(strcmp(s3->GetTxtDescriptor().c_str(), "3GC")==0);

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
}



void SnpAlignedTest::ImportCompressedBinary() {
	//2233466880   -> 10000101 00100000 00000000 00000000
	uint32_t smallArray[1];
	smallArray[0] = 2233466880;

	SnpAligned *snp1 = pool->GetSnp(4);
	snp1->SetIndividualCount(6);
	snp1->ImportCompressedBinSnp(0, 1, smallArray, 1, 2);
	BitSetType &g0=snp1->GetGenotype(1);
	//cout<<"\nBinary (g0): "<<g0<<"\n";
	CPPUNIT_ASSERT(g0[2]);
	CPPUNIT_ASSERT(!g0[0]);
	CPPUNIT_ASSERT(g0[3]);

	BitSetType &g1=snp1->GetGenotype(2);
	//cout<<"Binary (g1): "<<g1<<"\n";
	CPPUNIT_ASSERT(g1[0]);
	CPPUNIT_ASSERT(g1[5]);
	CPPUNIT_ASSERT(!g1[2]);
	CPPUNIT_ASSERT(!g1[3]);
	CPPUNIT_ASSERT(!g1[4]);

	BitSetType &g2=snp1->GetGenotype(3);
	//cout<<"Binary (g2): "<<g2<<"\n";
	CPPUNIT_ASSERT(!g2[0]);
	CPPUNIT_ASSERT(!g2[5]);
	pool->ReleaseSnp(snp1);
}	
	

void SnpAlignedTest::TestImporting() {
	
	SnpAligned *snp1 = pool->GetSnp(4);

	char *pico="312213\0";
	snp1->ImportSnp(0, 1, pico);
	BitSetType &g0=snp1->GetGenotype(1);
	//cout<<"\nTextBased (g0): "<<g0<<"\n";
	CPPUNIT_ASSERT(g0[1]);
	CPPUNIT_ASSERT(!g0[0]);
	CPPUNIT_ASSERT(g0[4]);

	BitSetType &g1=snp1->GetGenotype(2);
	//cout<<"TextBased (g1): "<<g1<<"\n";
	CPPUNIT_ASSERT(g1[2]);
	CPPUNIT_ASSERT(g1[3]);
	CPPUNIT_ASSERT(!g1[0]);
	CPPUNIT_ASSERT(!g1[1]);
	CPPUNIT_ASSERT(!g1[5]);

	BitSetType &g2=snp1->GetGenotype(3);
	//cout<<"TextBased (g2): "<<g2<<"\n";
	CPPUNIT_ASSERT(g2[0]);
	CPPUNIT_ASSERT(g2[5]);
	pool->ReleaseSnp(snp1);
}

void SnpAlignedTest::TestPunett() {
	SnpAligned *snp1 = pool->GetSnp(4);
	SnpAligned *snp2 = pool->GetSnp(4);

	char *map1 = "1211";			///<AA Aa AA AA
	char *map2 = "1332";			///<AA aa aa Aa
	
	snp1->ImportSnp(0, 1, map1);
	snp2->ImportSnp(1, 2, map2);
	
	
	SnpAligned *snp3 = snp1->punett(snp2);

	//We lost several genotypes because they are empty
	CPPUNIT_ASSERT(snp3->CountGenotypes() == 5);
	
	CPPUNIT_ASSERT(snp3->CountIndividuals(1) == 1);
	BitSetType &g0=snp3->GetGenotype(1);
	CPPUNIT_ASSERT(g0[0]);
	CPPUNIT_ASSERT(g0.count() == 1);
	
	
	BitSetType &g1=snp3->GetGenotype(2);
	CPPUNIT_ASSERT(g1[3]);
	CPPUNIT_ASSERT(g1.count() ==1);
	
	BitSetType &g2=snp3->GetGenotype(3);
	CPPUNIT_ASSERT(g2[2]);
	CPPUNIT_ASSERT(g2.count() == 1);
		

	BitSetType &g3=snp3->GetGenotype(4);
	CPPUNIT_ASSERT(g3[1]);
	CPPUNIT_ASSERT(g3.count() == 1);

	pool->ReleaseSnp(snp1);
	pool->ReleaseSnp(snp2);
	pool->ReleaseSnp(snp3);
}
void SnpAlignedTest::TestGenoLookup() {

	/* <UL><LI>000 - M
	 * <LI>001 - A
	 * <LI>010 - C
 	 * <LI>011 - G
	 * <LI>100 - T
	 */
	GenoLookup g1;
	g1.ParseBinaryGenotypes(33);			// 100001 - TA
	CPPUNIT_ASSERT(strcmp(g1.GetHeterozygote(), "TA")==0);
	CPPUNIT_ASSERT(strcmp(g1.GetHomozygote1(), "TT")==0);
	CPPUNIT_ASSERT(strcmp(g1.GetHomozygote2(), "AA")==0);

	g1.ParseBinaryGenotypes(26);			// 011010 - GC
	CPPUNIT_ASSERT(strcmp(g1.GetHeterozygote(), "GC")==0);
	CPPUNIT_ASSERT(strcmp(g1.GetHomozygote1(), "GG")==0);
	CPPUNIT_ASSERT(strcmp(g1.GetHomozygote2(), "CC")==0);
}
	
	
void SnpAlignedTest::TestUtil() {
	//cout<<"Testing basic usage patterns\n";
	SnpAligned *snps = pool->GetSnp(4);
	snps->SetIndividualCount(0);
	
	BitSetType &snp1 = snps->GetGenotype(1);
	snp1.append(5);
	CPPUNIT_ASSERT(snp1[0]);
	CPPUNIT_ASSERT(!snp1[1]);
	CPPUNIT_ASSERT(snp1[2]);
	BitSetType &snp2 = snps->GetGenotype(1);
	CPPUNIT_ASSERT(snp1 == snp2);
	BitSetType &snp3 = snps->GetGenotype(3);
	CPPUNIT_ASSERT(snp2 != snp3);
	CPPUNIT_ASSERT(snps->CountIndividuals(1) == 2);
	CPPUNIT_ASSERT(snps->CountIndividuals(3) == 0);

	pool->ReleaseSnp(snps);

	


}

void SnpAlignedTest::TestBitSetAssignment() {
	//First, I need to test the assignment of a bitset
	BitSetType b1;
	BitSetType b2;

	b1.resize(10);
	b1.set(1);
	b1.set(3);
	b1.set(5);
	b1.set(7);
	b1.set(9);
	b2=~b1;
	
	CPPUNIT_ASSERT(b2[0]);
	CPPUNIT_ASSERT(!b2[1]);
	CPPUNIT_ASSERT(b2[2]);
	CPPUNIT_ASSERT(!b2[3]);
	
}

void SnpAlignedTest::TestAssignment() {

	uint32_t smallArray[1];
	smallArray[0] = 2233466880;

	SnpAligned *snp1 = pool->GetSnp(4);
	SnpAligned *snp2=pool->GetSnp(4);

	
	snp1->SetIndividualCount(6);
	snp1->ImportCompressedBinSnp(0, 1, smallArray, 1, 2);

	//*snp2=*snp1;
	
	CPPUNIT_ASSERT(snp2->CountGenotypes() == 4);
	//CPPUNIT_ASSERT(snp2->CountIndividuals() == 6);


	BitSetType &g0=snp2->GetGenotype(1);
	//cout<<"\nBinary (g0): "<<g0<<"\n";
	CPPUNIT_ASSERT(g0[1]);
	CPPUNIT_ASSERT(!g0[0]);
	CPPUNIT_ASSERT(g0[4]);

	BitSetType &g1=snp2->GetGenotype(2);
	//cout<<"Binary (g1): "<<g1<<"\n";
	CPPUNIT_ASSERT(g1[2]);
	CPPUNIT_ASSERT(g1[3]);
	CPPUNIT_ASSERT(!g1[0]);
	CPPUNIT_ASSERT(!g1[1]);
	CPPUNIT_ASSERT(!g1[5]);

	BitSetType &g2=snp2->GetGenotype(3);
	//cout<<"Binary (g2): "<<g2<<"\n";
	CPPUNIT_ASSERT(g2[0]);
	CPPUNIT_ASSERT(g2[5]);
	pool->ReleaseSnp(snp2);
	pool->ReleaseSnp(snp1);
}



void SnpAlignedTest::TestLabels() {
	char snp1d[]="312";	
	char snp2d[]="213";

	SnpAligned *snp1=pool->GetSnp(4);
	SnpAligned *snp2=pool->GetSnp(4);
	
	snp1->ImportSnp(0, 1, snp1d);
	snp2->ImportSnp(1, 2, snp2d);
	
	CPPUNIT_ASSERT(strcmp(snp1->GetLabel(), "1") == 0);
	CPPUNIT_ASSERT(strcmp(snp2->GetLabel(), "2") == 0);

	SnpAligned *snp3 = snp1->punett(snp2);
	CPPUNIT_ASSERT(strcmp(snp3->GetLabel(), "1x2") == 0);

	pool->ReleaseSnp(snp1);
	pool->ReleaseSnp(snp2);
	pool->ReleaseSnp(snp3);
}

void SnpAlignedTest::TestMerge() {
	char snp1d[]="312";	
	char snp2d[]="213";

	SnpAligned *snp1=pool->GetSnp(4);
	SnpAligned *snp2=pool->GetSnp(4);
	
	snp1->ImportSnp(0, 1, snp1d);
	snp2->ImportSnp(1, 2, snp2d);

	snp1->Merge(snp2);


	BitSetType &g0=snp1->GetGenotype(1);
	//cout<<"\nTextBased (g0): "<<g0<<"\n";
	CPPUNIT_ASSERT(g0[1]);
	CPPUNIT_ASSERT(!g0[0]);
	CPPUNIT_ASSERT(g0[4]);

	BitSetType &g1=snp1->GetGenotype(2);
	//cout<<"TextBased (g1): "<<g1<<"\n";
	CPPUNIT_ASSERT(g1[2]);
	CPPUNIT_ASSERT(g1[3]);
	CPPUNIT_ASSERT(!g1[0]);
	CPPUNIT_ASSERT(!g1[1]);
	CPPUNIT_ASSERT(!g1[5]);

	BitSetType &g2=snp1->GetGenotype(3);
	//cout<<"TextBased (g2): "<<g2<<"\n";
	CPPUNIT_ASSERT(g2[0]);
	CPPUNIT_ASSERT(g2[5]);
	pool->ReleaseSnp(snp1);
	pool->ReleaseSnp(snp2);

}

}
}
