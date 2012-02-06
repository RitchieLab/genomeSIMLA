//
// C++ Interface: familynode
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSFAMILYNODE_H
#define GENETICSFAMILYNODE_H

#include "familymember.h"
#include "snpaligned.h"
#include <vector>
#include <iomanip>
#include "sibship.h"
namespace Genetics {

namespace Evaluation {
	class FamilyEvaluation;
}

using namespace std;

/**
	@brief Represent a family unit from a ped formatted file. 
	

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class FamilyNode{
public:
	/**
	 * This struct is used for each "run" to set up status and gentoype information that is used for the repository
	 */
	struct GenoData {
		GenotypeData *data;
		bool status;
	};

    FamilyNode(const char *id, BasicLog *log);

    ~FamilyNode();

	/**
	 * Returns the family member if one has been instantiated. Otherwise, it will create one
	 */
	FamilyMember *GetEntry(const char *id);

	/**
	 * @brief Iterate over the members and perform this evaluation
	 */
	bool PerformEvaluation( Evaluation::FamilyEvaluation *eval);


	 /**
	  * @brief Returns the number of individuals found inside
	  */
	uint GetMemberCount();


	/**
	 * @brief Returns the number of actual sibships have been identified within the pedigree
	 * It is not guaranteed that all sibships will participate in analyses. This is the raw number
	 * of sibships identified by common parental association.
	 */
	uint CountSibships();


	/**
	 * @brief Builds up the DSP pairs based on parental association. 
	 * In reality, this doesn't build pairs, but associates children with their biological parents so
	 * that sibs will be properly grouped together. It also keeps track of the numbers of affected/unaffected
	 * individuals
	 * @param TotalIndividualCount Indicate the size of the status Mask. This should be the total number of individuals seen so far
	 * @return dsp count. This indicates the number of pairs that were created. 
	 */
	uint BuildDsps(uint TotalIndividualCount);

	/**
	 * Returns the ID associated with the family node
	 */
	string GetID();

	/**
	 * @brief Produce a sib from the parents or an affected triad. 
	 * The resulting member is returned
	 */
	FamilyMember *CreateVirtualSib(FamilyMember *member);

	/**
	 * @brief Set the deletion status for an entire pedigree. 
	 */
	void MarkForDeletion(bool doMark);

	/**
	 * @brief Checks for deletion status
	 */
	bool MarkForDeletion();

	/**
	 * @brief returns the overall status for the family
	 */
	CaseControlStatus GetStatus();

	/**
	 * @brief Returns the overall scrambled status for the family
	 */
	CaseControlStatus GetScrambledStatus();

	/**
	 * @brief Produce a report describing the family and it's contents
	 */
	void GenerateReport(ostream *os);

	/**
	 * @brief Have the family append the relavent genodata to the vector
	 */
	uint AddGenoData(vector<GenotypeData *> &data, CaseControlStatus pgData[], uint &pgPosition,  bool doScrambleStatus /* = false */);
	
	//void SetFoldHome(uint i) { foldHome=i; }

	/**
 	 * @brief Returns the number of DSPs associated with a given family
	 */
	uint GetDspCount() { return dspCount;  }

	static bool expandAllAffecteds;						///<Indicate that we want all affecteds to be expanded with virtual siblings 
	static bool expandTrios;							///<Indicate that we want to expand families with valid parental data and 
														///<one affected child to a usable sibship
	
protected:
	string id;											///<The individual's id (from ped file)
	map <string, FamilyMember *> entries;				///<children associated with this member
	FamilyMemberPool *pool;								///<Pool used to acquire members
	ExtendedSibship sibships;							///<The siblings that have the same parents
	CaseControlStatus status;							///<The set of bitvectors representing the families status.This is valid AFTER the production of DSPs
	bool statusSet;
	bool markForDeletion;								///<Mark the whole family to be ignored
	BasicLog *pedLog;									///<This is written too when problems are encountered during setup
	//uint foldHome;										///<This is used in setting up the status masks for the folds.
	uint dspCount;										///<observed number of DSPs


};

inline
void FamilyNode::MarkForDeletion(bool doDelete) {
	markForDeletion=doDelete;
}

inline
bool FamilyNode::MarkForDeletion() {
	return markForDeletion;
}


inline
uint FamilyNode::CountSibships() {
	return sibships.size();
}	

inline
uint FamilyNode::GetMemberCount() {
	return entries.size();
}

inline
string FamilyNode::GetID() {
	return id;
}


inline
FamilyNode::FamilyNode(const char *id, BasicLog *pedLog) : id(id), statusSet(false), markForDeletion(false), pedLog(pedLog) {
	pool = FamilyMemberPool::Instance();
}

inline
FamilyNode::~FamilyNode() {
	FamilyMemberPool::Release();

	ExtendedSibship::iterator end=sibships.end();
	ExtendedSibship::iterator itr=sibships.begin();

	for (; itr!=end; itr++) {
		delete itr->second;
	}
	
}


inline
CaseControlStatus FamilyNode::GetStatus() {
	return status;
}

inline
CaseControlStatus FamilyNode::GetScrambledStatus() {
	cout<<"FamilyNode::GetScrambledStatus is dead code. please remove calling routines\n";
	assert(0);

}

inline
FamilyMember *FamilyNode::GetEntry(const char *id) {
	FamilyMember *member = NULL;
	string indID=this->id + "x" + string(id);

	map<string, FamilyMember *>::iterator end=entries.end();
	map<string, FamilyMember *>::iterator itr=entries.find(indID.c_str());
	if (atoi(id) > 0) {
		if (itr==end)	{
			member=pool->GetFamilyMember(this->id.c_str(), id);
			entries[indID]=member;
		}
		else	
			member=itr->second;
		
	}
	return member;
}

}

#endif
