//
// C++ Interface: sibshipfoldproduction
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTSIBSHIPFOLDPRODUCTION_H
#define MDRPDTSIBSHIPFOLDPRODUCTION_H
#include <vector>
#include "familyrepository.h"
namespace MdrPDT {
/**
	@author Eric Torstenson
*/
class FoldProduction{



struct FoldBucket {
	int dspCount;
	std::vector<FamilyRecord> families;
	FoldBucket() : dspCount(0) { }
	
	int AddFamily(FamilyRecord& f) {
		families.push_back(f);
		dspCount+=f.dspCount;
		return dspCount;
	}

	void WriteFoldToArray(int fold, char *&data) {
		std::vector<FamilyRecord>::iterator itr = families.begin();
		std::vector<FamilyRecord>::iterator end = families.end();
		
		while (itr != end ) {
			memset((void*)(data+itr->offset), (char)fold, itr->dspCount);
			itr++;
		}
	}

	bool operator<(const FoldBucket& other) const {
		return dspCount < other.dspCount; 
	}
};

typedef std::vector<FoldProduction::FoldBucket> BucketArray;


public:
    FoldProduction(uint foldCount)  : foldCount(foldCount), perfectScore(0), foldScore(-1), solutionsObserved(0) { }
	~FoldProduction() { }

	/**
	 * @brief performs the production, filling out the data pointed to by data
	 */
	void BuildXVFolds(char *data, std::vector<FamilyRecord> famRecords, int totalDSPs, Utility::Random& rnd);
	bool BuildXVFolds(int fIdx, int bIdx, BucketArray folds, Utility::Random& rnd);

	void GenerateReport(std::ostream& os);

	int GetDSPCount(int foldIndex) { return bestArrangement[foldIndex].dspCount;}
protected:
	int foldCount;
	std::vector<FoldBucket> buckets;
	std::vector<FamilyRecord> families;
	BucketArray bestArrangement;
	int perfectScore;
	int foldScore;
	int solutionsObserved;
};

}

#endif
