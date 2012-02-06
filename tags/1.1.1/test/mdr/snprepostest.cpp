//
// C++ Implementation: snprepostest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snprepostest.h"
#include "utility/utility.h"
#include "genetics/snprepository.h"
#include "evalmaxdifference.h"
#include "genetics/gtfileparserbuffered.h"
#include "genetics/gtlineparser.h"

CPPUNIT_TEST_SUITE_REGISTRATION( MDR::Test::SnpReposTest );

namespace MDR {

namespace Test {



SnpReposTest::SnpReposTest()
{
	pool=SnpPool::Instance();
}


SnpReposTest::~SnpReposTest()
{
	pool->Release();
}


void SnpReposTest::CompareSnps(SnpAligned *lhs, SnpAligned *rhs) {
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


void SnpReposTest::TestMultiLoading() {
/*	SnpRepository repo;

	char snp1d[]="22131221322223132213122132222313";			///<Line 0
	char snp2d[]="31233322223211323123332222321132";			///<Line 1			01111001010110101010011011110110
	char snp3d[]="32222332111111113222233211111111";			///<Line 9

	SnpAligned *s1=pool->GetSnp(4);
	s1->ImportSnp(1, 1, snp1d);
	
	SnpAligned *s2=pool->GetSnp(4);
	s2->ImportSnp(2, 2, snp2d);
	
	SnpAligned *s3=pool->GetSnp(4);
	s3->ImportSnp(3, 3, snp3d);
	
	repo.InitRepository(32, 0);
//	repo.InitRepository(32, 20);
	CPPUNIT_ASSERT(repo.GetIndividualCount() == 32);

	StringArray sar;
	sar.push_back("ExampleData10SNPs16Inds.txt");
	sar.push_back("ExampleData10SNPs16Inds.txt");
	sar.push_back("ExampleData10SNPs16Inds.txt");
	sar.push_back("ExampleData10SNPs16Inds.txt");


	//First, we need to initialize the repository. Let it know which files we want to use
	//repo.InitRepository("ExampleData10SNPs16Inds.txt");
	//repo.InitRepository("ExampleData10SNPs16Inds.txt");
	GtLineParserBP *lineParser=new GtLineParserBP();
	GtFileParserBuffered fileParser( lineParser, sar, 2);
	repo.ParseInputFile( &fileParser );
	CPPUNIT_ASSERT(repo.GetSnpCount() == 18);						//There are two that are skipped, so it's 18

	//repo.ParseBasePairTextFile("ExampleData10SNPs16Inds.txt");
	//repo.ParseBasePairTextFile("ExampleData10SNPs16Inds.txt");
	//repo.PostImport();
	
	SnpAligned *snp = repo.GetSnp(0);
	assert(snp != NULL);
	
	CompareSnps(s1, snp);
	CPPUNIT_ASSERT(snp!=NULL);
	snp=repo.GetSnp(10);
	CompareSnps(s2, snp);
	CPPUNIT_ASSERT(snp!=NULL);
	snp=repo.GetSnp(17);
	CPPUNIT_ASSERT(snp!=NULL);
	CompareSnps(s3, snp);

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
	*/
}

void SnpReposTest::TestMaxDiff() {
	BitSetType mask(16, 255ul);			//0000000011111111
	SnpRepository repo;
	repo.InitRepository(16, 10);
	repo.ParseEseBinGenofile("ExampleData10SNPs16Inds.bin");
	
	SnpRepository newRepo;
	newRepo.InitRepository(16, 0);
	newRepo.SetGrowby(10);

	CaseControlStatus status(mask);
	EvalMaxDifference md(6, 1);
	md.AppendStatus(status);
	repo.Evaluate(&newRepo, &md);
	
	CPPUNIT_ASSERT(newRepo.GetSnpCount() == 2);
	SnpAligned *s1 = newRepo.GetSnp((uint)0);
//	cout<<"Snp(0): "<<s1->GetLabel()<<"\n";
	CPPUNIT_ASSERT(strcmp(s1->GetLabel(), "5") == 0);
	SnpAligned *s2 = newRepo.GetSnp((uint)1);
//	cout<<"Snp(1): "<<s2->GetLabel()<<"\n";
	CPPUNIT_ASSERT(strcmp(s2->GetLabel(), "10") ==0);
//	pool->ReleaseSnp(s1);
//	pool->ReleaseSnp(s2);
}

void SnpReposTest::TestSelection() {
	BitSetType mask(10, 769ul);			//1100000001
	SnpRepository repo;
	repo.InitRepository(16, 10);
	repo.ParseEseBinGenofile("ExampleData10SNPs16Inds.bin");
	
	IsActive isactive(mask);
	SnpRepository newRepo;
	newRepo.InitRepository(16, 0);
	newRepo.SetGrowby(10);
	repo.Evaluate(&newRepo, &isactive);

	CPPUNIT_ASSERT(newRepo.GetSnpCount() == 3);
	SnpAligned *s1=newRepo.GetSnp((uint)0);		//2213122132222313
	CPPUNIT_ASSERT(s1->CountIndividuals(1) == 4);
	CPPUNIT_ASSERT(s1->CountIndividuals(2) == 8);
	CPPUNIT_ASSERT(s1->CountIndividuals(3) == 4);
	
	s1=newRepo.GetSnp(2);
	CPPUNIT_ASSERT(s1->CountIndividuals(1) == 8);
}


void SnpReposTest::TestPerfect() {
	BitSetType mask(16, 255ul);			//0000000011111111
	SnpRepository repo;
	repo.InitRepository(16, 10);
	repo.ParseEseBinGenofile("ExampleData10SNPs16Inds.bin");
	
	CaseControlStatus status(mask);
	IsPerfect isperfect(status);
	SnpRepository newRepo;
	newRepo.InitRepository(16, 0);
	newRepo.SetGrowby(10);
	repo.Evaluate(&newRepo, &isperfect);
	
	CPPUNIT_ASSERT(newRepo.GetSnpCount() == 1);
	SnpAligned *s1=newRepo.GetSnp((uint)0);
	CPPUNIT_ASSERT(s1->IsPerfect(mask) == 16);
}




void SnpReposTest::TestBinLoading() {
	SnpRepository repo;

	char snp1d[]="2213122132222313";			///<Line 0
	char snp2d[]="3123332222321132";			///<Line 1			01111001010110101010011011110110
	char snp3d[]="3222233211111111";			///<Line 9

	SnpAligned *s1=pool->GetSnp(4);

	//cout<<"\nImporting snp1:                     "<<snp1d<<"\n";
	s1->ImportSnp(100, 100, snp1d);
	
	//cout<<"\nImporting snp2:                     "<<snp2d<<"\n";
	SnpAligned *s2=pool->GetSnp(4);
	s2->ImportSnp(101, 101, snp2d);
	
	//cout<<"\nImporting snp9:                     "<<snp3d<<"\n";
	SnpAligned *s3=pool->GetSnp(4);
	s3->ImportSnp(102, 102, snp3d);

	//cout<<"\n";
	repo.InitRepository(16, 10);
	CPPUNIT_ASSERT(repo.GetSnpCount() == 10);
	CPPUNIT_ASSERT(repo.GetIndividualCount() == 16);

	repo.ParseEseBinGenofile("ExampleData10SNPs16Inds.bin");
	SnpAligned *snp = repo.GetSnp((uint)0);
	//cout<<"\nL: "<<s1->asciiGenotypes<<"\tR: "<<snp->asciiGenotypes<<"\n";
	assert(snp != NULL);
	
	CompareSnps(s1, snp);
	snp=repo.GetSnp(1);
	//cout<<"L: "<<s2->asciiGenotypes<<"\tR: "<<snp->asciiGenotypes<<"\n";
	CompareSnps(s2, snp);
	snp=repo.GetSnp(9);
	//cout<<"L: "<<s1->asciiGenotypes<<"\tR: "<<snp->asciiGenotypes<<"\n";
	CompareSnps(s3, snp);

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);

}


void SnpReposTest::TestMDRLoading() {
	SnpRepository repo;

	char snp1d[]="2213122132222313";			///<Line 0
	//char snp2d[]="3123332222321132";			///<Line 1			01111001010110101010011011110110
	char snp2d[]="1321112222123312";
	//char snp3d[]="3222233211111111";			///<Line 9
	char snp3d[]="1222211233333333";
	

	SnpAligned *s1=pool->GetSnp(4);
	s1->ImportSnp(100, 100, snp1d);
	
	SnpAligned *s2=pool->GetSnp(4);
	s2->ImportSnp(101, 101, snp2d);
	
	SnpAligned *s3=pool->GetSnp(4);
	s3->ImportSnp(102, 102, snp3d);
	
	repo.InitRepository(16, 0);
	
	CPPUNIT_ASSERT(repo.GetIndividualCount() == 16);	
	//First, we need to initialize the repository. Let it know which files we want to use
	//repo.InitRepository("ExampleData10SNPs16Inds.txt");

	GtLineParserMDR *lineParser=new GtLineParserMDR(1, false);
	GtFileParserBuffered *parser=new GtFileParserBuffered(lineParser, "ExampleData10SNPs16Inds.mdr");
	repo.ParseInputFile(parser, NULL);

	SnpAligned *snp = repo.GetSnp((uint)0);

		
	CPPUNIT_ASSERT(repo.GetSnpCount() == 9);	
	//repo.PostImport();
	snp = repo.GetSnp((uint)0);
	//cout<<"\nLine #1 "<<snp->toString();
	assert(snp != NULL);
	
	CaseControlStatus st;
	parser->GetStatusMask(&st);
	CPPUNIT_ASSERT(st.affected.count() == 8);
	CPPUNIT_ASSERT(st.unaffected.count()== 8);
	
	CompareSnps(s1, snp);
	snp=repo.GetSnp((uint)1);
	cout<<"Line #2 "<<snp->toString();
	CompareSnps(s2, snp);
	snp=repo.GetSnp((uint)8);
	cout<<"Line #9 "<<snp->toString();
	CompareSnps(s3, snp);

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
	delete parser;

}

void SnpReposTest::TestLoading() {
	SnpRepository repo;

	char snp1d[]="2213122132222313";			///<Line 0
	char snp2d[]="3123332222321132";			///<Line 1			01111001010110101010011011110110
	char snp3d[]="3222233211111111";			///<Line 9

	SnpAligned *s1=pool->GetSnp(4);
	s1->ImportSnp(100, 100, snp1d);
	
	SnpAligned *s2=pool->GetSnp(4);
	s2->ImportSnp(101, 101, snp2d);
	
	SnpAligned *s3=pool->GetSnp(4);
	s3->ImportSnp(102, 102, snp3d);
	
	repo.InitRepository(16, 10);
	
	CPPUNIT_ASSERT(repo.GetSnpCount() == 10);	
	CPPUNIT_ASSERT(repo.GetIndividualCount() == 16);

	//First, we need to initialize the repository. Let it know which files we want to use
	//repo.InitRepository("ExampleData10SNPs16Inds.txt");
	GtLineParserBP *lineParser=new GtLineParserBP(8,8,4);
	GtFileParserBuffered *parser=new GtFileParserBuffered(lineParser, "ExampleData10SNPs16Inds.txt");
	repo.ParseInputFile(parser, NULL);
	//repo.PostImport();
	SnpAligned *snp = repo.GetSnp((uint)0);
	//cout<<"\nLine #1 "<<snp->toString();
	assert(snp != NULL);
	
	CompareSnps(s1, snp);
	snp=repo.GetSnp(1);
	//cout<<"Line #2 "<<snp->toString();
	CompareSnps(s2, snp);
	snp=repo.GetSnp(8);
	//cout<<"Line #9 "<<snp->toString();
	CompareSnps(s3, snp);

	pool->ReleaseSnp(s1);
	pool->ReleaseSnp(s2);
	pool->ReleaseSnp(s3);
	delete parser;

}
}

}
