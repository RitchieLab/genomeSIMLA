//
// C++ Interface: pdtfold
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEPDTFOLD_H
#define ESEPDTFOLD_H
#include <vector>
#include "genetics/familynode.h"
#include <math.h>

namespace MDR {

using namespace Genetics;
using namespace std;

/**
@brief Encapsulates a single fold for a given portion of a cross validation test. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PdtFold {
public:
	struct Accuracy {
		uint falsePos;
		uint falseNeg;
		uint truePos;
		uint trueNeg;
		uint total;

		void CalculateTotal() {
			total=falsePos+falseNeg+truePos+trueNeg;
		}		

		void Reset() {
			falsePos 	= 0;
			falseNeg 	= 0;
			truePos 	= 0;
			trueNeg 	= 0;
			total 		= 0;
		}			
		Accuracy() : falsePos(0), falseNeg(0), truePos(0), trueNeg(0), total(0) {}
	};
	struct AccuracyEval {
		AccuracyEval() : trainingT(0.0), testingT(0.0), trainingD2(0), trainingD(0), testingD2(0), testingD(0)  { }
		void Reset() {
			classification.Reset();
			prediction.Reset();
			trainingD2	= 0;
			trainingD 	= 0;
			trainingT   = 0.0;
			testingD2	= 0;
			testingD	= 0;
			testingT	= 0.0;

		} 

		void CalculateTotal() {
			classification.CalculateTotal();
			prediction.CalculateTotal();
		}

		void CalcTrainingT(int D, int D2) {
			trainingD = D;
			trainingD2 = D2;
			trainingT = (float)trainingD / sqrt((float)trainingD2);
		}
		void CalcTestingT(int D, int D2) {
			testingD = D;
			testingD2 = D2;
			testingT = (float)testingD / sqrt((float)testingD2);
		}

		float ClassificationError() {
			return 100.0*(float)(classification.falsePos+classification.falseNeg)/((float)classification.total);
		}	

		float PredictionError() {
			return 100.0*(float)(prediction.falsePos+prediction.falseNeg)/((float)prediction.total);
		}

		float trainingT;
		float testingT;
		int trainingD2;
		int trainingD;	
		int testingD2;
		int testingD;
		Accuracy classification;
		Accuracy prediction;
		uint total;
	};

    PdtFold();
    ~PdtFold();

	/**
	 * @brief Initialize the statuses
	 */
	void Init(vector<FamilyNode*> families, uint individualCount);
	
	uint getDSPCount();						///<Returns the number of DSPs associated with this fold
	void SetDspCount(uint dsps);			///<Sets the number of DSps associated with the fold
	uint GetFamilyCount();					///<Returns the number of families associated with this node
	CaseControlStatus *GetStatuses();		///<Returns array of status (size of familyCount)

	/**
	 * @brief Calculate the denominator part of the T statistic
	 */
	uint CalculateD2(SnpAligned *model);
	void CalculateAccuracy(SnpAligned *model);
	AccuracyEval results;					///<Used to record the results of a given fold			

	CaseControlStatus GetOverallStatus();
	CaseControlStatus GetTrainingSet();
	void SetOverallStatus(CaseControlStatus& overall);
	
	PdtFold &operator=(const PdtFold&other);
	void ResetCache();
	static bool CacheD2;

protected:
	uint dspCount;							///<Indicates how many DSPs are represented by this fold
	uint familyCount;						///<Used to determine the size of the familyStatus array

	/**
	 * @brief affected/unaffected status for all individuals in this fold
	 */
	CaseControlStatus overallStatus;		///<The status for each member in the fold
	CaseControlStatus trainingSet;			///<The individuals used to train the algorithm
	/**
	 * @brief Array of statuses for each family in the fold
	 */
	CaseControlStatus *familyStatus;		///<Array of family status objects
	
	map<string, uint> d2Cache;
	uint cacheUses;

};


inline
PdtFold::PdtFold() : familyStatus(NULL), cacheUses(0) { }

inline
PdtFold::~ PdtFold()  {
	if (familyStatus)
		delete[] familyStatus;
}

inline
CaseControlStatus *PdtFold::GetStatuses() {
	return familyStatus;
}

inline
CaseControlStatus PdtFold::GetOverallStatus() {
	return overallStatus;
}

inline
void PdtFold::ResetCache() {
	if (CacheD2)
		d2Cache.clear();
}

inline
CaseControlStatus PdtFold::GetTrainingSet() {
	return trainingSet;
}

inline
uint PdtFold::getDSPCount() {
	return dspCount;
}


inline
void PdtFold::SetDspCount(uint count) {
	dspCount = count;
}

inline
uint PdtFold::GetFamilyCount() {
	return familyCount;
}

}

#endif
