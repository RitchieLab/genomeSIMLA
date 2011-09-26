//
// C++ Implementation: peddatatest
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <iostream>
#include <fstream>
#include "peddatatest.h"
#include "gtlineparsermdrpdt.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Genetics::Test::PedDataTest );

namespace Genetics {

namespace Test {

using namespace std;

PedDataTest::PedDataTest()
{
	familyRepo=NULL;
}


PedDataTest::~PedDataTest()
{
}

void PedDataTest::setUp() {
	LineParser lp;

	if (familyRepo == NULL) {
		familyRepo=FamilyRepository::Instance();
		PedParser pparser(familyRepo);
		lp.Parse("pedsample1.ped", &pparser);
	}
}

void PedDataTest::tearDown() {
	if (familyRepo)	{
		FamilyRepository::Release();
	}
	familyRepo=NULL;
}

void PedDataTest::TestFamilyArrangements() {
	FamilyNode *family = familyRepo->GetNode("1");
	
	CPPUNIT_ASSERT(familyRepo->GetFamilyCount() == 10);

	CPPUNIT_ASSERT(family->GetMemberCount() == 3);
	CPPUNIT_ASSERT(strcmp(family->GetID().c_str(), "1") ==0);
	FamilyMember *child = family->GetEntry("3");
	FamilyMember *father = family->GetEntry("1");
	CPPUNIT_ASSERT(child->GetFather() == father);
	CPPUNIT_ASSERT(father->GetMother() == NULL);
	CPPUNIT_ASSERT(father->GetFather() == NULL);
	//Just a few tests that the mom is the correct entry
	FamilyMember *mother = child->GetMother();
	CPPUNIT_ASSERT(strcmp(mother->GetID(), "2") == 0);
	CPPUNIT_ASSERT(strcmp(mother->GetFamilyID(), "1")==0);

}

/**
 * Tests:
 *		- Make sure that the individuals are loaded correctly
 * 		- Make sure that the families are set up properly
 */
void PedDataTest::TestBasicInput() {

	FamilyNode *family = familyRepo->GetNode("1");
	
	FamilyMember *child = family->GetEntry("3");
	FamilyMember *father = family->GetEntry("1");
	CPPUNIT_ASSERT(father->GetStatus() == 1);
	CPPUNIT_ASSERT(father->GetGenotypeValue(0) == 0);							//Test 0 0
	CPPUNIT_ASSERT(strcmp(father->GetGenotypeLabel(0).c_str(), "0 0") ==0);		
	CPPUNIT_ASSERT(father->GetGenotypeValue(1) == 2);							//Test 1 2
	CPPUNIT_ASSERT(strcmp(father->GetGenotypeLabel(1).c_str(), "1 2") == 0);
	CPPUNIT_ASSERT(father->GetGenotypeValue(2) == 2);							//Test 2 1
	CPPUNIT_ASSERT(strcmp(father->GetGenotypeLabel(2).c_str(), "1 2") == 0);
	CPPUNIT_ASSERT(father->GetGenotypeValue(3) == 1);							//Test 1 1
	CPPUNIT_ASSERT(strcmp(father->GetGenotypeLabel(3).c_str(), "1 1") == 0);
	CPPUNIT_ASSERT(father->GetGenotypeValue(6) == 3);							//Test 2 2
	CPPUNIT_ASSERT(strcmp(father->GetGenotypeLabel(6).c_str(), "2 2") == 0);
	CPPUNIT_ASSERT(father->GetGenotypeValue(13) == 2);							//Test the last entry

	//Just a few tests that the mom is the correct entry
	FamilyMember *mother = child->GetMother();
	CPPUNIT_ASSERT(mother->GetStatus() == 1);
	CPPUNIT_ASSERT(mother->GetGenotypeValue(0) == 2);
	CPPUNIT_ASSERT(mother->GetGenotypeValue(6) == 1);
	CPPUNIT_ASSERT(mother->GetGenotypeValue(13) == 1);
	
	//Some tests to verify that the child is the correct one
	CPPUNIT_ASSERT(child->GetStatus() == 2);
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 2);
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 1);
	CPPUNIT_ASSERT(child->GetGenotypeValue(13) == 2);

	family = familyRepo->GetNode("2");
	CPPUNIT_ASSERT(family->GetMemberCount() == 4);
	child = family->GetEntry("3");
	CPPUNIT_ASSERT(child->GetStatus() == 2);
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 2);
	CPPUNIT_ASSERT(child->GetGenotypeValue(13) == 3);		
	
	child = family->GetEntry("4");
	CPPUNIT_ASSERT(strcmp(child->GetFamilyID(), "2")==0);
	CPPUNIT_ASSERT(child->GetStatus() == 1);
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 3);
	CPPUNIT_ASSERT(child->GetGenotypeValue(13) == 1);		
	

	/*	
	WriteFamilies writeFamilies(&output, false);
	familyRepo->PerformEvaluation(&writeFamilies);	
	*/
}

void PedDataTest::TestErrorChecking() {
	//Test for various bad entries
	FamilyNode *family = familyRepo->GetNode("117");
	FamilyMember *child = family->GetEntry("3");
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 0);				//0 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(3) == 0);				//1 0
	CPPUNIT_ASSERT(child->GetGenotypeValue(7) == 0);				//2 0
	CPPUNIT_ASSERT(child->GetGenotypeValue(5) == 0);				//0 2

}

void PedDataTest::TestTrioFix() {
	FixTrios trioFix;
	familyRepo->PerformEvaluation(&trioFix);

	FamilyNode *family = familyRepo->GetNode("1");
	FamilyMember *child = family->GetEntry("16");
	cout<<"\n0 0 x 1 2 -> ??  x 0 0\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 0);
	cout<<"\n1 2 x 1 1 -> 1 1 & 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 1);
	cout<<"\n1 2 x 1 2 -> 2 2 & 1 1\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 1);	
	cout<<"\n1 1 x 1 1 -> 1 1 & 1 1\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(3) == 1);
	cout<<"\n2 2 x 1 2 -> 1 2 & 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(6) == 2);
	
	cout<<"\n1 1 x 1 2 -> 1 1 & 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(7) == 2);
	cout<<"\n1 1 x 1 2 -> 1 2 & 1 1\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(8) == 1);
	cout<<"\n2 2 x 1 2 -> 2 2 & 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(9) == 2);
	cout<<"\n2 2 x 1 2 -> 1 2 & 2 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(10) == 3);
	cout<<"\n1 2 x 1 2 -> 1 2 & 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(11) == 2);
	cout<<"\n1 2 x 1 2 -> 1 1 & 2 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(12) == 2);


	



	
}

void PedDataTest::TestMendelianErrorClearingSibs() {
	//Once we have cleaned up the values, they won't respond the same for the removal- have to restart
	FamilyRepository::Release();
	LineParser lp;
	familyRepo=FamilyRepository::Instance();
	PedParser pparser(familyRepo);
	lp.Parse("pedsample1.ped", &pparser);

	//This should now be propagated over to other children in the same family
	FamilyNode *family = familyRepo->GetNode("4");
	FamilyMember *child = family->GetEntry("4");
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 1);				//this becomes zero because the sibling is impossible
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 1);				//Sibling is impossible

	FixGenotypeErrors findErrors2(2, 6, &cout);
	familyRepo->PerformEvaluation(&findErrors2);
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 0);				//this becomes zero because the sibling is impossible
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 0);				//Sibling is impossible

}

void PedDataTest::TestMendelianErrorRemoves() {
	//Once we have cleaned up the values, they won't respond the same for the removal- have to restart
	FamilyRepository::Release();
	LineParser lp;
	familyRepo=FamilyRepository::Instance();
	PedParser pparser(familyRepo);
	lp.Parse("pedsample1.ped", &pparser);

	//Time for testing for 
	FamilyNode *family = familyRepo->GetNode("15");
	FamilyMember *child = family->GetEntry("1");
	CPPUNIT_ASSERT(!child->MarkForDeletion());

	FixGenotypeErrors findErrors3(3, 3, &cout);
	familyRepo->PerformEvaluation(&findErrors3);

	CPPUNIT_ASSERT(child->MarkForDeletion());

}

//Let's plant all of the mendal problems in family 15
void PedDataTest::TestMedelianError() {
	//Time for testing for 
	FamilyNode *family = familyRepo->GetNode("15");
	FamilyMember *child = family->GetEntry("1");
	
	//Check on the original data- this should get changed!
	CPPUNIT_ASSERT(!child->MarkForDeletion());
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 3);				//2 2  ! 1 1 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 1);				//1 1  ! 2 2 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 2);				//1 2  ! 1 1 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(3) == 1);				//1 1  ! 1 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(4) == 3);				//2 2  ! 1 2 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(5) == 2);				//1 2  ! 2 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(6) == 1);				//1 1  ! 2 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(7) == 3);				//2 2  ! 1 1 x 1 2

	FixGenotypeErrors findErrors0(0, 6, &cout);
	familyRepo->PerformEvaluation(&findErrors0);
	
	//Level 1 doesn't actual change anything!
	cout<<"\n2 2  ! 1 1 x 1 1\t "<<child->GetGenotypeLabel(0)<<"  ! "<<child->GetFather()->GetGenotypeLabel(0)<<" x "<<child->GetMother()->GetGenotypeLabel(0)<<"\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 3);				//2 2  ! 1 1 x 1 1
	cout<<"\n1 1  ! 2 2 x 1 1\t "<<child->GetGenotypeLabel(1)<<"  ! "<<child->GetFather()->GetGenotypeLabel(1)<<" x "<<child->GetMother()->GetGenotypeLabel(1)<<"\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 1);				//1 1  ! 2 2 x 1 1
	cout<<"\n1 2  ! 1 1 x 1 1\t "<<child->GetGenotypeLabel(2)<<"  ! "<<child->GetFather()->GetGenotypeLabel(2)<<" x "<<child->GetMother()->GetGenotypeLabel(2)<<"\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 2);				//1 2  ! 1 1 x 1 1
	cout<<"\n1 1  ! 1 2 x 2 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(3) == 1);				//1 1  ! 1 2 x 2 2
	cout<<"\n2 2  ! 1 2 x 1 1\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(4) == 3);				//2 2  ! 1 2 x 1 1
	cout<<"\n1 2  ! 2 2 x 2 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(5) == 2);				//1 2  ! 2 2 x 2 2
	cout<<"\n1 1  ! 2 2 x 2 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(6) == 1);				//1 1  ! 2 2 x 2 2
	cout<<"\n2 2  ! 1 1 x 1 2\t";
	CPPUNIT_ASSERT(child->GetGenotypeValue(7) == 3);				//2 2  ! 1 1 x 1 2

	FixGenotypeErrors findErrors1(1, 6, &cout);
	familyRepo->PerformEvaluation(&findErrors1);

	CPPUNIT_ASSERT(child->GetGenotypeValue(0) == 0);				//2 2  ! 1 1 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(1) == 0);				//1 1  ! 2 2 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(2) == 0);				//1 2  ! 1 1 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(3) == 0);				//1 1  ! 1 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(4) == 0);				//2 2  ! 1 2 x 1 1
	CPPUNIT_ASSERT(child->GetGenotypeValue(5) == 0);				//1 2  ! 2 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(6) == 0);				//1 1  ! 2 2 x 2 2
	CPPUNIT_ASSERT(child->GetGenotypeValue(7) == 0);				//2 2  ! 1 1 x 1 2
	




}

void PedDataTest::TestGenotypeData() {
 	char rawdata[]="0 0 1 1 1 2 0 1 1 1 2 2 1 2 2 1 1 1 0 0 2 0 0 2 1 0";
	char gt[4][4]={"0 0", "1 1", "1 2", "2 2"};
	int idx[]=  {0  ,1  ,2  ,0  ,1  ,3  ,2  ,2  ,1  ,0  ,0  ,0  ,0 };
	uint pos = 0;
	char geno[4];

	FamilyMember::GenoLkup lkupTable;

	GenotypeData data(&lkupTable);
	for (int i=0; i<13; i++) {
		pos = i*4;
		strncpy(geno, rawdata+pos, 3);
		geno[3]='\0'; 
		data.SetGenotype(i, geno);
		int val=data.GetGenotypeIndex(i);
		CPPUNIT_ASSERT(val == idx[i]);
		CPPUNIT_ASSERT(strcmp(data.GetGenotypeLabel(i).c_str(), gt[val]) == 0); 
	}	
	
}


}

}
