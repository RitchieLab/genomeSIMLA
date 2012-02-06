//
// C++ Interface: sibship
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSSIBSHIP_H
#define GENETICSSIBSHIP_H

namespace Genetics {

struct Sibship;

typedef map<std::string, Sibship *> ExtendedSibship;

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
struct Sibship {
/** 
 * @brief This structure is used to keep track of sibships (children that are associated with a given set of parents)
 * One of the primary roles of this structure is to recognize how these children should be weighted
 */
	Sibship(FamilyMember *father, FamilyMember *mother, uint indCount);
	Sibship(const Sibship &other);

	CaseControlStatus GetStatus();
	BitSetType GetPresence();

	/**
	 * We need to create the various masked genetoypes to handle missing 
 	 * data properly. So, we will create a mask for each possible combination. Then, when the
	 * p-tests are created, we will be getting the correct set of masks according to whom it is
	 * being "paired" with. 
	 */
	void AdjustWeights();
	void GenerateReport(ostream *os);


	/**
	 * @brief Insert a new child into the group
	 * @param child The new child to be added in
	 */
	void InsertChild(FamilyMember *child);

	
	vector<FamilyMember *> children;		///<Vector of children
	FamilyMember *father;					///<The father
	FamilyMember *mother;					///<The mother of the group

	uint affectedCount;						///<The number of affected individuals in the parental group
	uint unaffectedCount;					///<The Number of unaffected individuals in the parental group	
	uint unknownCount;						///<The number of individuals of Unknown status in the parental group
	
};


inline
Sibship::Sibship(FamilyMember *father, FamilyMember *mother, uint indCount) : father(father), mother(mother), affectedCount(0), unaffectedCount(0) { }		


inline
Sibship::Sibship(const Sibship &other) : father(other.father), mother(other.mother), affectedCount(other.affectedCount), unaffectedCount(other.unaffectedCount) {		
	int count = other.children.size();
	for (int i=0; i<count; i++)
		children.push_back(other.children[i]);

}


inline
CaseControlStatus Sibship::GetStatus() {
	CaseControlStatus status;
	for (uint i=0; i<children.size(); i++) {
		FamilyMember *child=children[i];

		uint indexCount=child->GetWeight();
		bool isAffected = child->IsAffected();
		//Going backwards to avoid resizing each time
		for (uint i=indexCount; i>0; i--) 
			status.SetStatus(child->GetIndividualIdx(i-1), isAffected);
	}
	status.affectedWeight=unaffectedCount;
	status.unaffectedWeight=affectedCount;
	return status;
}		


inline
BitSetType Sibship::GetPresence() {
	BitSetType presence;
	for (uint i=0; i<children.size(); i++) {
		FamilyMember *child=children[i];

		uint indexCount=child->GetWeight();

		//Going backwards to avoid resizing each time
		for (uint i=indexCount; i>0; i--) 
			presence[child->GetIndividualIdx(i-1)]=true;
	}
	return presence;
}

inline
void Sibship::AdjustWeights() {		
	uint count = children.size();
	for (uint i=0; i<count; i++) {
		FamilyMember *child1=children[i];
		for (uint j=i+1; j<count; j++) {
			FamilyMember *child2=children[j];
			child1->BuildMaskedGT(child2);
			child2->BuildMaskedGT(child1);
		}
	}
}

inline
void Sibship::GenerateReport(ostream *os) {
	stringstream details;
	stringstream parents;
	parents<<father->GetID()<<"x"<<mother->GetID();

	details<<"[ ";
	uint count = children.size();
	for (uint i=0; i<count; i++) {
		FamilyMember *child=children[i];

		if (child->IsAffected())
			details<<"*";
		else if (child->IsUnknownStatus()) 
			details<<"?";
		details<<child->GetID()<<" ";
	}
	details<<"]";
	*os<<setw(6)<<right<<parents.str()<<" "<<setw(16)<<left<<details.str();

}

inline
void Sibship::InsertChild(FamilyMember *child) {
	children.push_back(child);
	if (child->IsAffected())
		affectedCount++;
	else if (child->IsUnknownStatus()) 
		unknownCount++;
	else
		unaffectedCount++;
}



}

#endif
