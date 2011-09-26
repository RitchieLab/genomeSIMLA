//
// C++ Implementation: testpdttstat
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "testpdttstat.h"
#include "genetics/familynode.h"

CPPUNIT_TEST_SUITE_REGISTRATION( MDR::Test::TestPdtTStat );

namespace MDR {

namespace Test {

using namespace Genetics;

TestPdtTStat::TestPdtTStat() : families(NULL), pdtParser(NULL), fileParser(NULL), snps(NULL), evalPDT(NULL)
{
}


TestPdtTStat::~TestPdtTStat()
{
}
/**
 * 


	PdtFold *folds=p->GetFolds();
	p->GetStatusMask(&stat);
	evalPDT->SetFolds(folds, crossValidationCount, pgStatus, p->GetPgStatArrayCount());
	
	//Setup the overall status
	evalPDT->SetOverallStatus(stat);
	//Setup the trainer (for x-fold validation later on)
	evalSuite->SetTrainer(evalPDT);
*/
void TestPdtTStat::SetupRepository(const char *filename, bool doCreateVSibs) {
	EvalBalancedAccuracyPDT::Verbose=true;
	if (families == NULL)
		families = FamilyRepository::Instance();

	families->Purge(); 
	pdtParser=new GtLineParserMdrPdt(2, false, families, 1, 4 );
	FamilyNode::expandTrios = doCreateVSibs;
	pdtParser->SetPedigreeStream(&cout);

	fileParser=new GtFileParserBuffered(pdtParser, filename);
		
	snps = new SnpRepository();
	snps->ParseInputFile( fileParser, NULL );
	
	//Get the overall status mask 
	fileParser->GetStatusMask(&status);
	//FoldType st=pdtParser->GetStatArray();
	PdtFold *folds = pdtParser->GetFolds();

	pgStatus = pdtParser->GetPgStatArray();

	evalPDT = new EvalBalancedAccuracyPDT( 1, 1, 0.0, false, NULL, true);
	evalPDT->SetOverallStatus(status);
	evalPDT->SetFolds(folds, 1, pgStatus, pdtParser->GetPgStatArrayCount());

}

void TestPdtTStat::CleanupRepository() {
	if (pdtParser)
		delete pdtParser;
	pdtParser= NULL;

	fileParser=NULL;
	if (snps)
		delete snps;
	snps=NULL;
	
	if (evalPDT)
		delete evalPDT;
	evalPDT=NULL;	
	
	SnpPool::Instance()->Purge();

}

/**
 * @brief tests a set of models for expected affected/unaffected counts for each genotype
 */
void TestPdtTStat::VerifyGenotypeCounts(const char *modelID, uint affectedCount, uint unaffectedCount, uint expValues[] ) {
	cout<<status.affected<<"\n";
	cout<<status.unaffected<<"\n";
	
	CPPUNIT_ASSERT_EQUAL(affectedCount, status.affected.count());
	CPPUNIT_ASSERT_EQUAL(unaffectedCount, status.unaffected.count());

	bool passed = true;

	//OK, let's test a few models for the correct values:
	SnpAligned *model = snps->GetSnp(modelID);			//Get the top 1-SNP model we got when running the original
	cout<<"Model: "<<model->GetLabel()<<"\n";
//	cout<<"Present individuals:\n";
//	cout<<model->GetGenotype(0);
	
	uint gtCount = model->CountGenotypes() - 1;

	//Let's count the number of genotype 1s
	for (uint i = 0; i <gtCount; i++) { 
		BitSetType gt = model->GetGenotype(i+1);
		uint affecteds = (gt & status.affected).count();
		uint unaffecteds = (gt & status.unaffected).count();
		passed = passed && (expValues[i*2] == affecteds || expValues[i*2 + 1] == unaffecteds) ;
		
		if (passed) {
			CPPUNIT_ASSERT_EQUAL(expValues[i*2], affecteds);
			CPPUNIT_ASSERT_EQUAL(expValues[i*2 + 1], unaffecteds);
		}
		else
			cout<<i<<")"<<expValues[i*2]<<"/"<<affecteds<<"  "<<expValues[i*2+1]<<"/"<<unaffecteds<<" "; 
	}

	if (!passed) {
		cout<<"\n";
		evalPDT->EvaluateModel(model, 0);
		evalPDT->EvaluateVerbose( model );			
	}
	CPPUNIT_ASSERT(passed);
	model->ReduceInstanceCount();
}
/**
 * @brief tests a set of models for expected affected/unaffected counts for each genotype
 */
void TestPdtTStat::VerifyTScores(string models[], float tValues[], uint modelCount ) {
	//OK, let's test a few models for the correct values:
	SnpAligned *model = snps->GetSnp(models[0].c_str());
	evalPDT->EvaluateVerbose( model );
	ModelStatistics stats = evalPDT->EvaluateVerbose(model);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(tValues[0], stats.GetAvgTesting(), 0.001);
	
	for (uint i=0; i<modelCount; i++) {
		model=snps->GetSnp(models[i].c_str());
		stats = evalPDT->EvaluateVerbose(model);
		//evalPDT->EvaluateModel(model, i);
		//CPPUNIT_ASSERT_DOUBLES_EQUAL(tValues[i], model->GetLastMdEval(), 0.001);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(tValues[i], stats.GetAvgTesting(), 0.001);
		model->ReduceInstanceCount();

	}
}


void TestPdtTStat::TestModelsEpiAAU1x76() {
	cout<<"\nTesting EpiAAU.Model1.76.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 106, 90, 180, 204, 114, 106 };
	//Acquire the input parser from the configuration
	SetupRepository("EpiAAU.Model1.76.ped");
	VerifyGenotypeCounts("1", 400, 400, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={ 1.5428, 1.0776, 0.89517, 0.3410, 1.281025 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 0, 18, 109, 48, 0, 38, 94, 48, 0, 102, 105, 38, 0, 26, 92, 52, 0, 30 };
	VerifyGenotypeCounts("5x10", 400, 400, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = { 1.8329, 2.1193, 10.3441, 2.71413,  1.5475,  2.0412 };
	VerifyTScores(modelIDs, tValues, 5);
	CleanupRepository();
}

void TestPdtTStat::TestModelsEpiAAU1x100() {
	cout<<"\nTesting EpiAAU.Model1.100.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 106, 130, 192, 182, 102, 88 };
	//Acquire the input parser from the configuration
	SetupRepository("EpiAAU.Model1.100.ped");
	VerifyGenotypeCounts("2", 400, 400, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={ 1.9005, 1.9728, 1.8974, 1.2418, 0.06441 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 0, 32, 102, 56, 0, 14, 99, 38, 0, 96, 94, 60, 0, 22, 105, 64, 0, 18 };
	VerifyGenotypeCounts("5x10", 400, 400, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = { 2.0853, 2.5234, 9.5394, 2.1142, 2.1556, 2.1978  };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
}

void TestPdtTStat::TestModels2PCx31() {
	cout<<"\nTesting DSPs.Model2PC.31.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 46, 48, 107, 116, 47, 36};
	//Acquire the input parser from the configuration
	SetupRepository("DSPs.Model2PC.31.ped");
	VerifyGenotypeCounts("6", 200, 200, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={ 0.3511, 0.6975, 0.8006, 1.1921, 0.5774 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 8,11,9,20,29,13,17,25,77,59,10,18,39,17,7,29,4,8 };
	VerifyGenotypeCounts("5x10", 200, 200, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = { 1.9878, 1.7556, 5.9029, 1.9640, 1.5492, 2.0175 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
}

void TestPdtTStat::TestModels2PCx12() {
	cout<<"\nTesting DSPs.Model2PC.12.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 52, 44, 106, 103, 42, 53 };
	//Acquire the input parser from the configuration
	SetupRepository("DSPs.Model2PC.12.ped");
	VerifyGenotypeCounts("5", 200, 200, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={ 0.7071, 1.0000, 0.9258, 0.1302, 1.7179 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 8, 7, 15, 21, 29 ,16, 10, 23, 84, 54, 12, 26, 25, 16, 13, 27, 4, 10 };
	VerifyGenotypeCounts("5x10", 200, 200, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = { 2.3570, 1.9126, 5.7487, 2.0647, 1.5428, 2.4736 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
}

void TestPdtTStat::TestModels1x40() {
	cout<<"\nTesting DSPs.Model1.40.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={ 49, 67, 95, 88, 56, 45 };
	//Acquire the input parser from the configuration
	SetupRepository("DSPs.Model1.40.ped");
	VerifyGenotypeCounts("7", 200, 200, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={1.6165, 1.2309, 0.6882, 1.3333, 1.6876 };
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={ 0, 11, 48, 23, 0, 12, 55, 34, 0, 53, 44, 27, 0, 7, 53, 23, 0, 10 };
	VerifyGenotypeCounts("5x10", 200, 200, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = {1.9868, 2.9636, 9.6436, 2.3938, 2.4054, 2.4054 };
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
}


void TestPdtTStat::TestModels1x35() {
	cout<<"\nTesting DSPs.Model1.35.ped for known values\n";
	cout<<"1 SNP models\n";
	uint gtCounts[]={56, 54, 92, 86, 52, 60};
	//Acquire the input parser from the configuration
	SetupRepository("DSPs.Model1.35.ped");
	VerifyGenotypeCounts("5", 200, 200, gtCounts);
	string models[]={"1", "2", "3", "4", "5"};
	float tValues1Snp[]={1.3724, 2.3333, 2.037, 0.707, 1.2060};
	VerifyTScores(models, tValues1Snp, 5);
	cout<<"2 SNP models\n";
	uint gtCounts2Snp[]={0,11, 56, 28, 0, 15, 36, 23, 0, 39, 56, 24, 0, 12, 52, 28, 0, 20};
	VerifyGenotypeCounts("5x10", 200, 200, gtCounts2Snp);
	string modelIDs[] = {"1x2", "1x10", "5x10", "7x8", "8x9", "9x10"};
	float tValues[] = {2.3635, 1.6733, 9.8489, 1.8856, 2.6499, 3.1497};
	VerifyTScores(modelIDs, tValues, 5);

	CleanupRepository();
	
}


}
}
