//
// C++ Implementation: sibshipfoldproduction
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "foldproduction.h"
#include <iomanip>
#include <sstream>
namespace MdrPDT {

using namespace std;

/**
 * @brief performs the production, filling out the data pointed to by data
 */
void FoldProduction::BuildXVFolds(char *data, std::vector<FamilyRecord> famRecords, int totalDSPs, Utility::Random& rnd) {
	families = famRecords;
	
	//Time to shuffle the families...kinda. We want them sorted according to
	//their DSP count, but we want them mixed up otherwise.So, we need to 
	//shuffle them first, then do the sort (which doesn't guarantee order
	//past the DSP Count)
	random_shuffle(families.begin(), families.end(), rnd);


	//Now just a sort, using our < operator
	sort(families.begin(), families.end());
	
	buckets.clear();
	buckets.resize(foldCount);

	//This should guarantee that the very first one is considered worth saving
	foldScore = totalDSPs;
	//A perfect score has exactly the same number in each bucket
	perfectScore = totalDSPs % foldCount;
	if (perfectScore > 0 && perfectScore < famRecords[famRecords.size() - 1].dspCount)
		perfectScore = famRecords[famRecords.size() - 1].dspCount;
	
	//Start the search
	BuildXVFolds(0, 0, buckets, rnd);


	//Populate the array
	std::vector<FoldBucket>::iterator itr = bestArrangement.begin();
	std::vector<FoldBucket>::iterator end = bestArrangement.end();

	int fold = 0;
	while (itr != end) {
		itr->WriteFoldToArray(fold++, data);
		itr++;
	}

}

bool FoldProduction::BuildXVFolds(int fIdx, int bIdx, BucketArray folds, Utility::Random& rnd) {
	bool perfectSolutionFound = false;
	uint attempts = 0;
	assert(fIdx < (int)families.size());
	

	FamilyRecord currRecord = families[fIdx++]; 
	folds[bIdx].AddFamily(currRecord);
	sort(folds.begin(), folds.end());

	//Check to see if we have exhausted our population; 
	if (fIdx < (int)families.size()) {
		while (!perfectSolutionFound && (int)attempts<foldCount) {
			//If we failed last time, there is no point in pushing the same number again
			if (attempts == 0 || !(folds[attempts].dspCount == folds[attempts-1].dspCount))
				//check to make sure we don't have a perfect solution 
				perfectSolutionFound = BuildXVFolds(fIdx, attempts, folds, rnd);
			attempts++;
		}
	}
	else {
		//We have reached a potential solution, so let's evaluate it and see if we should keep it
		//The vector should already be sorted in ascending order...so...
		int score = folds[foldCount - 1].dspCount - folds[0].dspCount;

//		for (uint i = 1; i< foldCount; i++) 
//			score += folds[i].dspCount - folds[i-1].dspCount;
		//Let's check to see if it's better than a previous score
		if (score < foldScore) {
			foldScore = score;
			bestArrangement = folds;
		}
		

		perfectSolutionFound = score <= perfectScore;
	}
	return perfectSolutionFound;		
}



void FoldProduction::GenerateReport(std::ostream& os) {
	os<<"\n\nFold Selection Details: \n";
	os<<setw(25)<<"Number of Folds: "<<foldCount<<"\n";
	os<<setw(25)<<"Perfect Score: "<<perfectScore<<"\n";
	os<<setw(25)<<"Final Score: "<<foldScore;
	if (foldScore == perfectScore)		
		os<<" Perfect!\n";
	else
		os<<"\n"; 
	os<<setw(25)<<"Solutions Observed: "<<solutionsObserved<<"\n";

	std::vector<FoldBucket>::iterator itr = bestArrangement.begin();
	std::vector<FoldBucket>::iterator end = bestArrangement.end();

	int foldID=0;
	while (itr != end) {
		stringstream ss;
		os<<"Fold "<<++foldID<<": Total Families: "<<itr->families.size()<<"\n";
		os<<"Individual Contributions: ";
		ss<<"Pedigree IDs:      ";
		for (uint n=0; n<itr->families.size(); n++) {
			ss<<" "<<itr->families[n].pedigreeID;
			os<<" "<<itr->families[n].dspCount;
		}
		os<<"\n"<<ss.str()<<"\n";
		itr++;
	}

	os<<endl;	
	
	//We used to display the individuals, but I'm not sure how we'd do that right now....
	//Maybe a series of X's where the individual was part of the slice?	
}

}
