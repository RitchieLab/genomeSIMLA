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


/**
 * The overall status tells us which individuals are going to be used for the whole search. 
 * We need to find everyone that is  in common with overall and the inverse of the testing set. 
 */
void PdtFold::SetOverallStatus(CaseControlStatus& overall) {
/*	cout<<"Setting Overall Status:\n";
	cout<<"    Previous Values: "<<trainingSet.total.count()
		<<" = "<<trainingSet.affected.count()<<"+"<<trainingSet.unaffected.count()<<"\n";
	cout<<"    Total: "<<setw(8)<<overallStatus.total.count()<<" "<<overallStatus.total<<"\n";
	cout<<"    Aff  : "<<setw(8)<<overallStatus.affected.count()<<" "<<overallStatus.affected<<"\n";
	cout<<"    Unaff: "<<setw(8)<<overallStatus.unaffected.count()<<" "<<overallStatus.unaffected<<"\n";
	*/	
	trainingSet = CaseControlStatus(overallStatus);
	trainingSet.Flip();
	trainingSet.total&=overall.total;
	trainingSet.affected&=overall.affected;
	trainingSet.unaffected&=overall.unaffected;	

/*	cout<<"    Training Status:\n";
	cout<<"    Total: "<<setw(8)<<trainingSet.total.count()<<" "<<trainingSet.total<<"\n";
	cout<<"    Aff  : "<<setw(8)<<trainingSet.affected.count()<<" "<<trainingSet.affected<<"\n";
	cout<<"    Unaff: "<<setw(8)<<trainingSet.unaffected.count()<<" "<<trainingSet.unaffected<<"\n";
	*/
	//cout<<"Fold Description: \n\t"<<trainingSet.affected<<"\n\t"<<trainingSet.unaffected<<"\n";
		
}
void PdtFold::CalculateAccuracy(SnpAligned *model) {
//void PdtFold::CalculateAccuracy(uint hrCells[], uint cellCount, SnpAligned *model) {
	uint clAffCount		= 0;		///<Count of affected individuals
	uint clUnaffCount 	= 0;		///<Count of unaffected individuals
	uint prAffCount     = 0;
	uint prUnaffCount	= 0;

	uint gtCount = model->CountGenotypes();
	results.classification.Reset();
	results.prediction.Reset();

	BitSetType hrCells = model->GetHrCells();

	//cout<<"X Training Set: "<<trainingSet.total<<"\n";
	//cout<<"X Testing Set:  "<<overallStatus.total<<"\n";
	//cout<<"X Testing Aff:  "<<overallStatus.affected<<"\n";


	BitSetType curGenotype;
	//no need to worry about families here
	for (uint i=1; i<gtCount; i++) {
		curGenotype = model->GetGenotype(i);
		prAffCount   = (curGenotype & overallStatus.affected).count();
		prUnaffCount = (curGenotype & overallStatus.unaffected).count();
		clAffCount   = (curGenotype & trainingSet.affected).count();
		clUnaffCount = (curGenotype & trainingSet.unaffected).count();

		if (hrCells[i]) {
		//Check to see if it's high risk
		//if (cellIdx < cellCount && i == hrCells[cellIdx]) {
			results.prediction.falsePos+=prUnaffCount;
			results.prediction.truePos+=prAffCount;
			results.classification.falsePos+=clUnaffCount;
			results.classification.truePos+=clAffCount;
			//cellIdx++;
		} else {
			results.prediction.falseNeg+=prAffCount;
			results.prediction.trueNeg+=prUnaffCount;
			results.classification.falseNeg+=clAffCount;
			results.classification.trueNeg+=clUnaffCount;
		}
			
	}
	results.CalculateTotal();

}

uint PdtFold::CalculateD2(SnpAligned *model) {
	uint affCount		= 0;		///<Count of affected individuals
	uint unaffCount 	= 0;		///<Count of unaffected individuals
	uint sumD2 			= 0;		///<The sum of D squared
	int localD;						///<The local D for a given family

	string label;

	BitSetType hrCells = model->GetHrCells();

	if (CacheD2) {
		stringstream lbl;
		lbl<<hrCells;
/*		lbl<<model->GetLabel();
		for (uint i =0; i<cellCount; i++) 
			lbl<<hrCells[i];
	*/
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
		size_t i 				= 0;		///<Index used in iterating over genotypes
		localD = 0;
		BitSetType curGenotype;
		for (i = hrCells.find_first(); i != BitSetType::npos; i=hrCells.find_next(i)) {
		//for (i=0; i<cellCount; i++) {
			curGenotype=model->GetGenotype(i);
			
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
