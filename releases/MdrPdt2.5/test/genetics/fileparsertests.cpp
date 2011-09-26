//
// C++ Implementation: fileparsertests
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fileparsertests.h"
#include "genetics/gtfileparserbuffered.h"
#include "genetics/gtlineparser.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Genetics::Test::FileParserTests );

namespace Genetics {

namespace Test {

using namespace Genetics::Parser;

FileParserTests::FileParserTests()
{
	pool=SnpPool::Instance();
}


FileParserTests::~FileParserTests()
{
	SnpPool::Release();
}


void FileParserTests::CompareSnps(SnpAligned *lhs, SnpAligned *rhs) {
	cout<<"\nlhs: "<<lhs->toString()<<"\nrhs: "<<rhs->toString()<<"\n";
	BitSetType &l=lhs->GetGenotype(1);
	BitSetType &r=rhs->GetGenotype(1);
	cout<<"\nl: "<<l<<"\nr: "<<r<<"\n";
	CPPUNIT_ASSERT((l&r) == l);

	l=lhs->GetGenotype(2);
	r=rhs->GetGenotype(2);
	cout<<"l: "<<l<<"\nr: "<<r<<"\n";
	CPPUNIT_ASSERT((l&r) == l);
	
	l=lhs->GetGenotype(3);
	r=rhs->GetGenotype(3);
	cout<<"l: "<<l<<"\nr: "<<r<<"\n\n\n";
	CPPUNIT_ASSERT((l&r) == l);
}



void FileParserTests::TestMdrLineParser() {
	char line[]="1 0 1 2 0 1 2 0 1 2 0";
	char data1[]="1111";
	char data2[]="2222";
	char data3[]="3333";

	
	GtLineParserMDR lineParser(1, false);

	lineParser.ParseLine(line, 0);
	lineParser.ParseLine(line, 1);
	lineParser.ParseLine(line, 2);
	lineParser.ParseLine(line, 3);

	SnpAligned *s1=pool->GetSnp(4);
	SnpAligned *s2=pool->GetSnp(4);
	SnpAligned *s3=pool->GetSnp(4);
	SnpAligned *s4=pool->GetSnp(4);
	s1->ImportSnp(10, 10, data1);
	s2->ImportSnp(11, 11, data2);
	s3->ImportSnp(12, 12, data3);
	s4->ImportSnp(13, 13, data1);

	SnpAligned *snp=lineParser.GetSnp(0);
	CompareSnps(snp, s1);
	pool->ReleaseSnp( snp );
	
	snp=lineParser.GetSnp(1);
	CompareSnps(snp, s2);
	pool->ReleaseSnp( snp );

	snp=lineParser.GetSnp(2);
	CompareSnps(snp, s3);
	pool->ReleaseSnp( snp );

	snp=lineParser.GetSnp(3);
	CompareSnps(snp, s4);
	pool->ReleaseSnp( snp );

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
	pool->ReleaseSnp(s4);
	
}

void FileParserTests::TestMdrLoading() {
	char snp1d[]="2213122132222313";			///<Line 0
	char snp2d[]="1321112222123312";			///<Line 1			01111001010110101010011011110110
	char snp3d[]="2132213323212221";
	char snp4d[]="2123333233332312";
	char snp5d[]="1223122122323322";
	char snp6d[]="1212333232223113";
	char snp7d[]="2222222222222222";
	char snp8d[]="2123132133122212";
	char snp9d[]="2123132122213123";
	char snp10d[]="1222211233333333";


	SnpAligned *s1=pool->GetSnp(4, 16);
	s1->ImportSnp(100, 100, snp1d);
	
	SnpAligned *s2=pool->GetSnp(4, 16);
	s2->ImportSnp(101, 101, snp2d);
	
	SnpAligned *s3=pool->GetSnp(4, 16);
	s3->ImportSnp(102, 102, snp3d);

	SnpAligned *s4=pool->GetSnp(4, 16);
	s4->ImportSnp(103, 103, snp4d);

	SnpAligned *s5=pool->GetSnp(4, 16);
	s5->ImportSnp(104, 104, snp5d);

	SnpAligned *s6=pool->GetSnp(4, 16);
	s6->ImportSnp(105, 105, snp6d);

	SnpAligned *s7=pool->GetSnp(4, 16);
	s7->ImportSnp(106, 106, snp7d);

	SnpAligned *s8=pool->GetSnp(4, 16);
	s8->ImportSnp(107, 107, snp8d);

	SnpAligned *s9=pool->GetSnp(4, 16);
	s9->ImportSnp(108, 108, snp9d);

	SnpAligned *s10=pool->GetSnp(4, 16);
	s10->ImportSnp(104, 104, snp10d);


	GtLineParserMDR *mdrParser=new GtLineParserMDR(1, false);
	GtFileParserBuffered parser(mdrParser, "ExampleData10SNPs16Inds.mdr");
	parser.Open();
	SnpAligned *snp;

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s1);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();	
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s2);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s3);
	pool->ReleaseSnp( snp );
		
	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s4);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s5);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s6);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s7);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s8);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s9);
	pool->ReleaseSnp( snp );

	snp=parser.NextSnp();
	CPPUNIT_ASSERT(snp);
	CompareSnps(snp, s10);
	pool->ReleaseSnp( snp );


	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
	pool->ReleaseSnp(s4);
	pool->ReleaseSnp(s5);
	pool->ReleaseSnp(s6);
	pool->ReleaseSnp(s7);
	pool->ReleaseSnp(s8);
	pool->ReleaseSnp(s9);
	pool->ReleaseSnp(s10);


}

}

}
