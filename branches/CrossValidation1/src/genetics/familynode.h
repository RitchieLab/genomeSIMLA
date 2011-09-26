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
	uint CountSibships();
	/**
	 * @brief Builds up the DSP pairs based on parental association. 
	 * In reality, this doesn't build pairs, but associates children with their biological parents so
	 * that sibs will be properly grouped together. It also keeps track of the numbers of affected/unaffected
	 * individuals
	 * @param TotalIndividualCount Indicate the size of the status Mask. This should be the total number of individuals seen so far
	 * @param createVirtSibs allow the user to control the production of virtual siblings
	 * @return dsp count. This indicates the number of pairs that were created. 
	 */
	uint BuildDsps(uint TotalIndividualCount, bool createVirtSibs);

	/**
	 * Returns the ID associated with the family node
	 */
	string GetID();

	/**
	 * @brief Produce a sib from the parents or an affected triad. 
	 * The resulting member is returned
	 */
	FamilyMember *CreateVirtualSib(FamilyMember *member);

	void MarkForDeletion(bool doMark);
	bool MarkForDeletion();

	/**
	 * @brief returns the overall status for the family
	 */
	CaseControlStatus *GetStatus();

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
	uint AddGenoData(vector<GenotypeData *> &data, CaseControlStatus ccData[], CaseControlStatus pgData[], uint &ccPosition, uint &pgPosition,  bool doScrambleStatus /* = false */);
protected:
	string id;											///<The individual's id (from ped file)
	map <string, FamilyMember *> entries;				///<children associated with this member
	FamilyMemberPool *pool;								///<Pool used to acquire members
	ExtendedSibship sibships;							///<The siblings that have the same parents
	CaseControlStatus status;							///<The set of bitvectors representing the families status.This is valid AFTER the production of DSPs
	bool statusSet;
	bool markForDeletion;								///<Mark the whole family to be ignored
	BasicLog *pedLog;									///<This is written too when problems are encountered during setup
	

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
CaseControlStatus *FamilyNode::GetStatus() {
	stringstream os;
	
	if (markForDeletion)
		return NULL;
	if (statusSet) {
		return &status;
	}
	statusSet=true;
	//Otherwise, we need to put it all together
	map<string, FamilyMember *>::iterator end=entries.end();
	map<string, FamilyMember *>::iterator itr=entries.begin();

	os<<"Setting up status for family: "<<GetID()<<"\n";
	if (itr==end)
		os<<"No entries found\n";
	
	for (; itr!=end; itr++) {
		FamilyMember *child=itr->second;
		child->AppendStatus(status);
	}
	os<<"\n\t"<<status.affected.count()<<"\t"<<status.unaffected.count()<<"\t"<<status.total.count()<<"\n";
	if (pedLog)
		pedLog->Write( os.str().c_str() );
	return &status;
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
