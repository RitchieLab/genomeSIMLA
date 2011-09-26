//
// C++ Interface: sibshipfoldproduction
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESESIBSHIPFOLDPRODUCTION_H
#define ESESIBSHIPFOLDPRODUCTION_H
#include "pdtfold.h"
#include "genetics/familyrepository.h"
#include <vector>
namespace MDR {

using namespace std;
using namespace Genetics;



/**
@brief Put together folds for x-validation

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SibshipFoldProduction{
public:
struct FamilyRecord {
	FamilyNode *family; 
	uint dspCount;

	FamilyRecord() : family(NULL), dspCount(0) { }
	FamilyRecord(FamilyNode *family) : family(family) {
		dspCount = family->GetDspCount();
	}
	
	//we want these to be sorted descending
	bool operator< (const FamilyRecord& other) const {
		return dspCount > other.dspCount;
	}
};

struct FoldBucket {
	CaseControlStatus * fold;	
	uint dspCount;
	vector<FamilyNode *> families;

	FoldBucket() : fold(NULL), dspCount(0) { }
	FoldBucket(CaseControlStatus * fold) : fold(fold), dspCount(0) { }
	BitSetType GetMask();

	uint AddFamily(FamilyRecord &f) {
		families.push_back(f.family);
		dspCount+=f.family->GetDspCount();
		return dspCount;
	}	


	//We want this one to be sorted in ascending order
	bool operator< (const FoldBucket& other) const {
		return dspCount < other.dspCount;
	}
};

	typedef vector<SibshipFoldProduction::FoldBucket> BucketArray;

    SibshipFoldProduction(uint foldCount, FamilyRepository *families, uint individualCount);

    ~SibshipFoldProduction();

	/**
	 * @brief Builds the cross validation folds. Returns a pointer to the array of folds. 
	 * The production object destroys them when it goes away. 
	 */
	PdtFold* BuildXVFolds(uint individualCount);

	void GenerateReport(ostream *os);
	CaseControlStatus GetOverallStatus();
protected:
	/**
	 * @brief The recusrive search component.
	 * This should almost exhaust all possible combinations if it can't find a perfect solution.
	 * Right now, a perfect score is 0. However, there will be some sets where the perfect score
	 * is actually higher than that. A score represents the difference from one bucket to the next. 
	 */
	bool BuildXVFolds(uint famPos, uint bucket, BucketArray folds);

	FamilyRepository *families;				///<Repository of family data
	uint foldCount;							///<The number of folds we are splitting everything into
	uint foldScore; 						///<This is used to short circuit the search. 
	PdtFold *folds;				///<The array of folds. This is built up once we find a solution (though it won't always be the optimal solution until the end)
	BucketArray bestArrangement;
	vector<FamilyRecord> familyData; 
	uint perfectScore;						///<This is an estimate of the perfect score
	uint solutionsObserved;					///<Used to record the number of solutions were considered before the current solution was found
	CaseControlStatus overallStatus;

};

inline
SibshipFoldProduction::SibshipFoldProduction(uint foldCount, FamilyRepository *families, uint individualCount) : families(families), foldCount(foldCount), solutionsObserved(0) {
	folds = new PdtFold[foldCount];
	foldScore = 0;
	foldScore--;
	//Yeah, this sort of defeats the typedef	
}

inline
SibshipFoldProduction::~SibshipFoldProduction() {
//	delete[] folds;
}

inline
CaseControlStatus SibshipFoldProduction::GetOverallStatus() {
	return overallStatus;
}

}

#endif
