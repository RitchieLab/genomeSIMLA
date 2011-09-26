//
// C++ Interface: fixtrios
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFIXTRIOS_H
#define GENETICS_EVALUATIONFIXTRIOS_H
#include "familyrepoevaluationbasic.h"
#include "utility/types.h"
namespace Genetics {

namespace Evaluation {

#define UNAFFECTED 1

/**
@brief Evaluate members of a family for being trios and perform the necessary work to make them ready for PDT analysis
 */
class FixTrios : public FamilyRepoEvaluationBasic {
public:
	FixTrios();
	~FixTrios();

	bool operator()(FamilyNode *family, uint position);
	/**
	 * @brief Performs an evaluation to determine the usefulness of a given family member
	 */
	bool operator()(FamilyMember *individual, uint position);
	void Report(ostream *os);	
protected:
	uint totalTriosSeen;
	uint invalidTrios;
	uint validTrios;
	uint parentsRemoved;
	uint childrenCreated;
	uint nonTrioFamiliesEncountered; 
};

inline
bool FixTrios::operator()(FamilyNode *family, uint position) {
	this->family = family;
	repoCount++;
	totalTriosSeen++;
	if (family->GetMemberCount() == 3) {
		cout<<"Evaluating Family: "<<family->GetID()<<"\tMember Count: "<<family->GetMemberCount()<<"\t";
		return family->PerformEvaluation(this);
	}
	nonTrioFamiliesEncountered++;
	return false;
}

inline
void FixTrios::Report(ostream *os) {
	*os<<"Trio Evaluation Summary:";
	*os<<"\nTotal Trios Observed:       "<<totalTriosSeen;
	*os<<"\nInvalid Trios Encountered:  "<<invalidTrios;
	*os<<"\nValid Trios:                "<<validTrios;
	*os<<"\nParents Removed:            "<<parentsRemoved;
	*os<<"\nChildren Created:           "<<childrenCreated;
	*os<<"\nnonTrioFamiliesEncountered: "<<nonTrioFamiliesEncountered;
	*os<<"\n";
}

class ReportFamilies : public FamilyRepoEvaluationBasic {
public:
	ReportFamilies(ostream *os) : os(os) { }
	~ReportFamilies() {}

	bool operator()(FamilyMember *member, uint position);
protected:
	ostream *os;
};

/**
@brief Writes family data in a given family repository to a file

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class WriteFamilies : public FamilyRepoEvaluationBasic {
public:
	/**
	 * @brief write the output to os and return the value, returnVal
	 */
	WriteFamilies(ostream *os, bool returnVal) : rval(returnVal),os(os) {}

	WriteFamilies() {}
	~WriteFamilies(){};

	/**	
  	 * @brief Determine if the locus will get written to the stream
	 * @note By default, all loci are written
	 */
	virtual bool VerifyLocus(uint locus){ return true; };	

	/**
	 * @brief For each family member, we are going to write their data to the stream (along with the meaningless stuff)
	 */
	bool Evaluate(FamilyMember *family, uint position);
protected:
	bool rval;
	ostream *os;
};

class WriteSelectedLoci : public WriteFamilies {
public:
	WriteSelectedLoci(ostream *os, BitSetType &loci, bool returnVal);
	WriteSelectedLoci(ostream *os, SnpAligned *model, bool returnVal);
	~WriteSelectedLoci();
	/**	
  	 * @brief Determine if the locus will get written to the stream
	 * @note By default, all loci are written
	 */
	virtual bool VerifyLocus(uint locus);	
protected:
	bool SetupModel(SnpAligned *model);

	BitSetType validLoci;
};

inline
FixTrios::FixTrios() : totalTriosSeen(0), invalidTrios(0), validTrios(0), parentsRemoved(0), childrenCreated(0), nonTrioFamiliesEncountered(0) {}


inline
FixTrios::~FixTrios() {}

inline
WriteSelectedLoci::WriteSelectedLoci(ostream *os, BitSetType &loci, bool returnVal) 
	: WriteFamilies(os, returnVal), validLoci(loci) { 
}

inline
WriteSelectedLoci::WriteSelectedLoci(ostream *os, SnpAligned *model, bool returnVal) 
	: WriteFamilies(os, returnVal) {
		SetupModel(model);
}

inline
WriteSelectedLoci::~WriteSelectedLoci() { }





		


}

}

#endif
