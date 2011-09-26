//
// C++ Implementation: testpdttstatreal
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "testpdttstatreal.h"
CPPUNIT_TEST_SUITE_REGISTRATION( ESE::Test::TestPdtTStatReal );

namespace ESE {

namespace Test {

TestPdtTStatReal::TestPdtTStatReal()
{
}


TestPdtTStatReal::~TestPdtTStatReal()
{
}



void TestPdtTStatReal::TestModelsSCZ() {
	cout<<"\nTesting SCZ_D5.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 224, 298, 382, 319, 91, 80 };
	//Acquire the input parser from the configuration
	SetupRepository("SCZ_d5.ped", false);
	VerifyGenotypeCounts("31", 758, 758, gtCounts);
	string models[]={"31", "48", "29", "12", "2" };
	float tValues1Snp[]={ 3.673, 3.201, 3.105, 0.9461, 0.5467 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 75, 45, 73, 95, 64, 45, 87, 76, 193, 147, 56, 105, 31, 62, 56, 72, 46, 34 };
	VerifyGenotypeCounts("60x71", 758, 758, gtCounts2Snp);
	string modelIDs[] = { "60x71", "38x42", "40x47", "28x37", "37x42" };
	float tValues[] = {4.6790, 3.12863, 2.1527, 2.4874, 2.1980 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();

}

void TestPdtTStatReal::TestModelsALZ() {
	cout<<"\nTesting Fam_alldsps.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={373, 705, 777, 623, 263, 85};
	//Acquire the input parser from the configuration
	SetupRepository("Fam_alldsps.ped");
	VerifyGenotypeCounts("31", 1482, 1482, gtCounts);
	string models[]={"1", "10", "22", "31", "48"};
	float tValues1Snp[]={ 0.8969,  1.2971, 3.3191, 8.1782, 2.07688 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 12, 21, 19, 38, 22, 6, 80, 182, 187, 150, 65, 17, 221, 371, 428, 326, 129, 52 };
	VerifyGenotypeCounts("21x31", 1482, 1482, gtCounts2Snp);
	string modelIDs[] = { "21x31", "38x42", "40x47", "28x37", "37x42" };
	float tValues[] = {8.21401, 1.9718, 2.6383, 1.8209, 2.4361 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();

}



void TestPdtTStatReal::TestModelsBPADCIT() {
	cout<<"\nTesting BPADCIT_DISC_GRINtrios.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={76, 64, 4, 27, 30, 19};
	//Acquire the input parser from the configuration
	SetupRepository("BPADCIT_DISC_GRINtrios.ped");
	VerifyGenotypeCounts("55", 113, 113, gtCounts);
	string models[]={"1", "15", "25", "40", "55", "65"};
	float tValues1Snp[]={0.3922, 1.1339, 2.7217, 1.2185, 4.2710, 1.27920};
	VerifyTScores(models, tValues1Snp, 6);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 54,35,3,21,28,19,21,20,0,5,2,0,0,8,1,1 };
	VerifyGenotypeCounts("51x55", 113, 113, gtCounts2Snp);
	string modelIDs[] = {"6x15", "21x31", "13x17", "35x40", "51x55" };
	float tValues[] = {1.1209, 1.5000, 1.5690, 1.5213, 5.096369 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
	
}


}

}
