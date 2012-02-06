//
// C++ Implementation: pdtfold
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pdtfold.h"

namespace MDR {

void PdtFold::Init(vector<FamilyNode *>families, uint individualCount) {
	familyCount = families.size();
	if (familyStatus)
		delete[] familyStatus;

	//Start with an empty status
	overallStatus.Reset();
	familyStatus = new CaseControlStatus[familyCount];
	

	//Iterate over each member and capture their status
	for (uint i=0; i< familyCount; i++) {
		familyStatus[i] = families[i]->GetStatus();
		familyStatus[i].Resize(individualCount);
		overallStatus.AppendStatus(familyStatus[i]);
	}
}	

void PdtFold::SetOverallStatus(CaseControlStatus& overall) {
	trainingSet = CaseControlStatus(overallStatus);
	trainingSet.Flip();
	trainingSet.affected&=overall.affected;
	trainingSet.unaffected&=overall.unaffected;	
		
}

void PdtFold::CalculateAccuracy(uint hrCells[], uint cellCount, SnpAligned *model) {
	uint affCount		= 0;		///<Count of affected individuals
	uint unaffCount 	= 0;		///<Count of unaffected individuals

	uint gtCount = model->CountGenotypes();
	uint cellIdx = 0;
	
	results.falsePos = 0;
	results.falseNeg = 0;
	results.truePos  = 0;
	results.trueNeg  = 0;
	BitSetType curGenotype;
	//no need to worry about families here
	for (uint i=1; i<gtCount; i++) {
		curGenotype = model->GetGenotype(i);
		affCount = (curGenotype & overallStatus.affected).count();
		unaffCount = (curGenotype & overallStatus.unaffected).count();

		//Check to see if it's high risk
		if (cellIdx < cellCount && i == hrCells[cellIdx]) {
			results.falsePos+=unaffCount;
			results.truePos+=affCount;
			cellIdx++;
		} else {
			results.falseNeg+=affCount;
			results.trueNeg+=unaffCount;
		}
			
	}
	results.CalculateTotal();

}

uint PdtFold::CalculateD2(uint hrCells[], uint cellCount, SnpAligned *model) {
	uint i 				= 0;		///<Index used in iterating over genotypes
	uint affCount		= 0;		///<Count of affected individuals
	uint unaffCount 	= 0;		///<Count of unaffected individuals
	uint sumD2 			= 0;		///<The sum of D squared
	int localD;						///<The local D for a given family

	CaseControlStatus *lastFold = &familyStatus[familyCount];

	//Walk through the family masks
	for (CaseControlStatus *currCC=familyStatus; currCC<lastFold; currCC++) {
		localD = 0;
		BitSetType curGenotype;
		for (i=0; i<cellCount; i++) {
			curGenotype=model->GetGenotype(hrCells[i]);
			
			affCount = (currCC->affected & curGenotype).count();
			unaffCount = (currCC->unaffected & curGenotype). count();
			localD+=affCount - unaffCount;
		}
		sumD2+=(localD*localD);

	}					
	
	return sumD2;
}

}
