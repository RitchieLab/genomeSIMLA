//
// C++ Implementation: sibshipfoldproduction
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "sibshipfoldproduction.h"

namespace MDR {


void SibshipFoldProduction::GenerateReport(ostream *os) {
	*os<<"\n\nFold Selection Details: \n";
	*os<<setw(25)<<"Number of Folds: "<<foldCount<<"\n";
	*os<<setw(25)<<"Perfect Score: "<<perfectScore<<"\n";
	*os<<setw(25)<<"Final Score: "<<foldScore;
	if (foldScore == perfectScore)		
		*os<<" Perfect!\n";
	else
		*os<<"\n"; 
	*os<<setw(25)<<"Solutions Observed: "<<solutionsObserved<<"\n";

	for (uint i=0; i<foldCount; i++) {
		stringstream ss;
		*os<<"Fold "<<i+1<<": Total Families: "<<bestArrangement[i].families.size()<<"\n";
		*os<<"DSP Contributions: ";
		 ss<<"Pedigree IDs:      ";
		for (uint n=0; n<bestArrangement[i].families.size(); n++) {
			ss<<" "<<bestArrangement[i].families[n]->GetID();
			*os<<" "<<bestArrangement[i].families[n]->GetDspCount();
		}
		*os<<"\n"<<ss.str()<<"\n";
	}
	*os<<"\nIndividuals: "<<endl;	
	
	for (uint i=0; i<foldCount; i++) 
		*os<<"Fold:       "<<i+1<<" "<<folds[i].GetOverallStatus().total<<endl;

}


//Gotta be careful, folds will be copied every single call
bool SibshipFoldProduction::BuildXVFolds( uint fIdx, uint bIdx, BucketArray folds) {
	bool perfectSolutionFound = false;
	uint attempts = 0;
	assert(fIdx < familyData.size());

	FamilyRecord currRecord = familyData[fIdx++]; 
	folds[bIdx].AddFamily(currRecord);
	sort(folds.begin(), folds.end());

	//Check to see if we have exhausted our population; 
	if (fIdx < familyData.size()) {
		while (!perfectSolutionFound && attempts<foldCount) {
			//If we failed last time, there is no point in pushing the same number again
			if (attempts == 0 || !(folds[attempts].dspCount == folds[attempts-1].dspCount))
				//check to make sure we don't have a perfect solution 
				perfectSolutionFound = BuildXVFolds(fIdx, attempts, folds);
			attempts++;
		}
	}
	else {
		//We have reached a potential solution, so let's evaluate it and see if we should keep it
		//The vector should already be sorted in ascending order...so...
		uint score = folds[foldCount - 1].dspCount - folds[0].dspCount;

		for (uint i = 1; i< foldCount; i++) 
			score += folds[i].dspCount - folds[i-1].dspCount;
		//Let's check to see if it's better than a previous score

		for (uint i=0; i<foldCount; i++) {
			stringstream ss;
		}

		if (score < foldScore) {
			foldScore = score;
			bestArrangement = folds;
		}


		perfectSolutionFound = score <= perfectScore;
	}
	return perfectSolutionFound;		
}

PdtFold *SibshipFoldProduction::BuildXVFolds(uint individualCount) {
	FamilyRepository::Iterator itr = families->GetIterator();
	FamilyNode *family = itr.GetNext();
	uint overallDSPs = 0;

	familyData.clear();
	//Build up the family vector 
	while (family) {
		if (!family->MarkForDeletion() && family->GetDspCount() > 0) {
			overallDSPs+= family->GetDspCount();
			FamilyRecord rec(family);	
			familyData.push_back(rec);
		}
		family = itr.GetNext();

	}
	
	//OK, we are going to shuffle the contents so that they are mixed up before we 
	//perform the "sort". Keep in mind the sort doesn't guarantee any order when the DSP counts
	//are the same, so we should be able to rely on a different order going into the sort should 
	//yield different orders
	random_shuffle(familyData.begin(), familyData.end(), Utility::Random::globalGenerator);

	//The record should be setup so that it sorts in descending order
	sort(familyData.begin(), familyData.end());
	BucketArray buckets(foldCount);

	perfectScore = overallDSPs % foldCount * 2;
	if (perfectScore > 0 && perfectScore < (2 * familyData[familyData.size() - 1].dspCount ) )
		perfectScore = 2 * familyData[familyData.size() - 1].dspCount;

	//Perform the search. 
	bool perfectSolution = BuildXVFolds(0, 0, buckets);
	overallStatus.Reset();
	if (!perfectSolution) {
		//We might need to do something special to make sure we are happy with the score
	}
	for (uint i=0; i<foldCount; i++)  {
		folds[i].Init(bestArrangement[i].families, individualCount);
		CaseControlStatus s = folds[i].GetOverallStatus();
		overallStatus.AppendStatus(s);
		folds[i].SetDspCount(bestArrangement[i].dspCount);
	}

/*	cout<<"    Total: "<<setw(8)<<overallStatus.total.count()<<" "<<overallStatus.total<<"\n";
	cout<<"    Aff  : "<<setw(8)<<overallStatus.affected.count()<<" "<<overallStatus.affected<<"\n";
/*	cout<<"    Unaff: "<<setw(8)<<overallStatus.unaffected.count()<<" "<<overallStatus.unaffected<<"\n";
*/
	for (uint i=0; i<foldCount; i++)  {

/*		cout<<"    Total: "<<setw(8)<<overallStatus.total.count()<<" "<<overallStatus.total<<"\n";
		cout<<"     +Aff : "<<setw(8)<<overallStatus.affected.count()<<" "<<overallStatus.affected<<"\n";
/*		cout<<"    Unaff: "<<setw(8)<<overallStatus.unaffected.count()<<" "<<overallStatus.unaffected<<"\n";
*/
		folds[i].SetOverallStatus(overallStatus);
	}

	
	return folds;
}


}
