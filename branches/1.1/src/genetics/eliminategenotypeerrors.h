//
// C++ Interface: fixeliminategenotypeerrors
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFIXELIMINATEGENOTYPEERRORS_H
#define GENETICS_EVALUATIONFIXELIMINATEGENOTYPEERRORS_H
#include "familyrepoevaluationbasic.h"
#include <iostream>
#include <vector>
namespace Genetics {

namespace Evaluation {

using namespace std;

/**
@brief Evaluates the genotypes for a given family to ensure that the genotypes of the children are valid and therefore reliable for consideration in the statistic.

Any affected children that don't match their parents are rendered as unknown.

If configuration is set to do this, the genotypes for all family members will be set to unknown if one of the affected children 
<table border="0"><tr><th bgcolor="#dadada">Value </th><th bgcolor="#dadada"> Action Taken</th></tr>
<tr><td bgcolor="#ffffff" class="twikiFirstCol"> 0 </td><td bgcolor="#ffffff"> Report only </td></tr>
<tr><td bgcolor="#eaeaea" class="twikiFirstCol"> 1 </td><td bgcolor="#eaeaea"> Report and zero out that child's genotype information at that loci </td></tr>

<tr><td bgcolor="#ffffff" class="twikiFirstCol"> 2 </td><td bgcolor="#ffffff"> Report, zero out the values at that loci for all members of that family </td></tr>
<tr><td bgcolor="#eaeaea" class="twikiFirstCol"> 3 </td><td bgcolor="#eaeaea"> If a single child exceeds a certain threshold, the whole family will be removed from the data set - regardless, all members will be zero'd out for the missing loci) </td></tr>
</table>
<H2> Reported information </h2>
<table border="0"><tr><th bgcolor="#dadada" class="twikiFirstCol">Entry Frequency</th><th bgcolor="#dadada"> Entry Details</th></tr>

<tr><td bgcolor="#ffffff" class="twikiFirstCol"> Per Error Encountered </td><td bgcolor="#ffffff"> Family ID, Parent IDs, Child ID, All three individual's genotype values </td></tr>
<tr><td bgcolor="#eaeaea" class="twikiFirstCol"> Per Child Evaluated </td><td bgcolor="#eaeaea"> Summary of each child's bad and missing loci </td></tr>
<tr><td bgcolor="#ffffff" class="twikiFirstCol"> Per Family </td><td bgcolor="#ffffff"> Summary of total number of loci across the whole family (if loci 1 is missing in 3 members, it will be counted only once)** </td></tr>
<tr><td bgcolor="#eaeaea" colspan="2" class="twikiFirstCol">  </td></tr>

<tr><td bgcolor="#ffffff" class="twikiFirstCol">  </td><td bgcolor="#ffffff"> <em>*What exactly is a family? Is it anyone with the same family ID or only those with the same parents?</em> </td></tr>
</table>

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class FixGenotypeErrors : public FamilyRepoEvaluationBasic {
public:
	/**
	 * @brief Construction
	 * @param level Indicates behavior when we encounter errors
	 * @param threshold Used to indicate the point at which a family is deemed invalid
	 */
    FixGenotypeErrors(uint level, uint threshold, ostream *os);
	FixGenotypeErrors();

    ~FixGenotypeErrors();


	bool Evaluate(FamilyNode *family, uint position);
	/**
	 * @brief Identifies any invalid genotypes for each child (with parent information) and zeros them out
	 */
	bool Evaluate(FamilyMember *individual, uint position);
	 
	 /**
	  * Generate complete sumary
	  */
	void Report(ostream *os);	

	/**
	 * @brief used to reset the affected loci vector
	 */
	void Reset();

	/**
	 * @brief Array used to keep up with which loci need to be purged, if we exceed the threshold
	 */
	vector<uint> affectedCounts;
protected:
	/**
	 * This is used internally to resize the number of loci are associated with the family
	 */
	void SetLociCount(uint number);
	string HandleLoci(uint loci, FamilyMember *member, uint father, uint mother, uint child);
	ostream *report;
	uint badData;
	uint totalPedErrors;
	uint overallErrors;								///<The total error count for all families seen
	uint level;
	uint threshold;
	BitSetType affectedLoci;
	uint maxErrorCount;
	uint totalGenotypes;							///<Indicates how many genotypes observed (max for a family)
	

};


class ClearBadFamilies : public FamilyEvaluation {
public:
	ClearBadFamilies() {}
	~ClearBadFamilies() {}

	/**
	 * @brief Identifies any invalid genotypes for each child (with parent information) and zeros them out
	 */
	bool Evaluate(FamilyMember *individual, uint position);

	 /**
	  * Generate complete sumary
	  */
	void Report(ostream *os);	
};


class StripBadLoci : public FamilyRepoEvaluationBasic {
public:
	StripBadLoci(vector<uint> *affCounts, uint mode, uint threshold) : mode(mode), affCounts(affCounts), threshold(threshold), report(affCounts->size(), 0) {	}

	StripBadLoci() : mode(0), affCounts(NULL), threshold(0) {	}

	~StripBadLoci() {}

	/**
	 * @brief Identifies any invalid genotypes for each child (with parent information) and zeros them out
	 */
	bool Evaluate(FamilyMember *individual, uint position);

	/**
	 * @brief Make sure the details are set up properly
	 */
	void Init();
	 
	 /**
	  * @brief Generate complete sumary
	  */
	void Report(ostream *os);	

	/**
	 * @brief used to reset the affected loci vector
	 */
	void Reset();

protected:
	uint mode;								///<What mode are we running in (3 means act. Otherwise, just report)
	vector<uint> *affCounts;				///<This is used to determine how many errors have occured for each locus
	uint threshold;							///<This is used to determine if the error count is bad enough to clear
	vector<uint> report;					///<Used to tally the number actually stripped
};


inline
bool ClearBadFamilies::Evaluate(FamilyMember *individual, uint position) {
	individual->MarkForDeletion(true);
	return false;
}

inline
void ClearBadFamilies::Report(ostream *os) {
	
}

class StripSibGT : public FamilyEvaluation {
public:
	StripSibGT(BitSetType *validLoci) : affectedLoci(validLoci) {}
	~StripSibGT() {}
	
	/**
	 * @brief Identifies any invalid genotypes for each child (with parent information) and zeros them out
	 */
	bool Evaluate(FamilyMember *individual, uint position);	

	 /**
	  * Generate complete sumary
	  */
	void Report(ostream *os);	
protected:
	BitSetType *affectedLoci;
};

inline
FixGenotypeErrors::FixGenotypeErrors() : report(NULL), totalPedErrors(0), level(1), threshold(75) {
	badData = 0;
	totalPedErrors = 0;
	maxErrorCount = 0;
	totalGenotypes = 0;							///<Indicates how many genotypes observed (max for a family)
	overallErrors=0;
}

inline
FixGenotypeErrors::FixGenotypeErrors(uint level, uint threshold, ostream *os) : report(os), totalPedErrors(0), level(level), threshold(threshold) {
	badData = 0;
	totalPedErrors = 0;
	maxErrorCount = 0;
	totalGenotypes = 0;							///<Indicates how many genotypes observed (max for a family)
	overallErrors=0;
}

inline
FixGenotypeErrors::~FixGenotypeErrors() {}



}

}

#endif
