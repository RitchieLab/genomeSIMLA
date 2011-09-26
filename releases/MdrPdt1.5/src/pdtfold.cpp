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

bool PdtFold::CacheD2 = false;



/**
 * Gets data from the other object...except for the cache. The cache is assumed to be reset
 */
PdtFold &PdtFold::operator=(const PdtFold& other) {
	dspCount=other.dspCount;
	familyCount=other.familyCount;
	overallStatus=other.overallStatus;
	trainingSet=other.trainingSet;
	familyStatus = new CaseControlStatus[familyCount];
	for (uint i=0; i<familyCount; i++) 
		familyStatus[i]=other.familyStatus[i];
	ResetCache();
	return *this;
}
		

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
	uint clAffCount		= 0;		///<Count of affected individuals
	uint clUnaffCount 	= 0;		///<Count of unaffected individuals
	uint prAffCount     = 0;
	uint prUnaffCount	= 0;

	uint gtCount = model->CountGenotypes();
	uint cellIdx = 0;
	results.classification.Reset();
	results.prediction.Reset();


	BitSetType curGenotype;
	//no need to worry about families here
	for (uint i=1; i<gtCount; i++) {
		curGenotype = model->GetGenotype(i);
		prAffCount   = (curGenotype & overallStatus.affected).count();
		prUnaffCount = (curGenotype & overallStatus.unaffected).count();
		clAffCount   = (curGenotype & trainingSet.affected).count();
		clUnaffCount = (curGenotype & trainingSet.unaffected).count();

		//Check to see if it's high risk
		if (cellIdx < cellCount && i == hrCells[cellIdx]) {
			results.prediction.falsePos+=prUnaffCount;
			results.prediction.truePos+=prAffCount;
			results.classification.falsePos+=clUnaffCount;
			results.classification.truePos+=clAffCount;
			cellIdx++;
		} else {
			results.prediction.falseNeg+=prAffCount;
			results.prediction.trueNeg+=prUnaffCount;
			results.classification.falseNeg+=clAffCount;
			results.classification.trueNeg+=clUnaffCount;
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

	string label;

	if (CacheD2) {
		stringstream lbl;
		lbl<<model->GetLabel();
		for (uint i =0; i<cellCount; i++) 
			lbl<<hrCells[i];
	
		label = lbl.str();
	
		if (d2Cache.size() > 0) {
			map<string, uint>::iterator itr=d2Cache.find(label);
			if (itr != d2Cache.end())   {
				cacheUses++;
				return itr->second;	
			}
		}
	}

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

	if (CacheD2)		
		d2Cache[label]=sumD2;
	return sumD2;
}

}
