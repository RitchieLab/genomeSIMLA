//
// C++ Implementation: familynode
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "familynode.h"
#include "familyevaluation.h"
#include <iomanip>

namespace Genetics {
using namespace Evaluation;
using namespace std;

bool FamilyNode::expandAllAffecteds = false;
bool FamilyNode::expandTrios        = false;


bool FamilyNode::PerformEvaluation( FamilyEvaluation *eval) {
	map <string, FamilyMember *>::iterator itr = entries.begin();
	map <string, FamilyMember *>::iterator end = entries.end();

	FamilyMember *member;
	uint count=0;
	for (; itr != end; itr++) {
		member=itr->second;

		//Delete the whole family if the return tells you so
		if (eval->Evaluate(member, count++)) {
			//delete member;
			//entries.erase(itr);
		}
	}
	return entries.size() < 1;
}

void FamilyNode::GenerateReport(ostream *os) {
	ExtendedSibship::iterator itr = sibships.begin();
	ExtendedSibship::iterator end = sibships.end();	
	stringstream stbuff;
	if (os) {
		uint affectedCount = 0;
		uint unaffectedCount = 0;
		uint balancedCount = 0;

		for (; itr!=end; itr++)  {
			affectedCount+=itr->second->affectedCount;
			unaffectedCount+=itr->second->unaffectedCount;
			balancedCount+=((itr->second->affectedCount * itr->second->unaffectedCount)*2);
			itr->second->GenerateReport(&stbuff);
		}
		*os<<setw(12)<<right<<GetID();

		if (markForDeletion)
			*os<<setw(8)<<0<<setw(8)<<affectedCount<<setw(8)<<unaffectedCount;
		else
			*os<<setw(8)<<balancedCount<<setw(8)<<affectedCount<<setw(8)<<unaffectedCount;

		*os<<" : "<<stbuff.str()<<"\n";
	}
}


/**
@brief Removes trio parents from the family and creates the virtual control.
Basically, if the family member is a parent with 1 child, they will be removed from the family.
If the member is a child of a trio group, create the virtual control and add it to the family.
New child will have an id of 25 higher than the original affected child.
*/ 
FamilyMember *FamilyNode::CreateVirtualSib(FamilyMember *member) {
	ostream *os = NULL;
	if (pedLog)
		os = pedLog->GetStream();

	FamilyMember *father = member->GetFather();
	FamilyMember *mother = member->GetMother();
	FamilyMember *newChild = NULL;
	
	//If both parents aren't missing, then we can't do anything with them
	if ((father == NULL && mother == NULL) || member->MarkForDeletion() )  {
		return false;
	}
	
	if (member->IsAffected()) {
		//create the copy 
		char id[16];
		sprintf(id, "%lu",  entries.size()+15);
		
		GenotypeData *childData 	= member->GetGenotypeData();
		GenotypeData *fatherData 	= father->GetGenotypeData();
		GenotypeData *motherData	= mother->GetGenotypeData();

		if (fatherData->CountGenotypesPresent() * motherData->CountGenotypesPresent() == 0) {
			if (os)
				*os<<"Sibship: "<<member->GetFamilyID()<<":"<<father->GetID()<<"x"<<mother->GetID()<<"\t\tUnable to create virtual sibling due to lack of parental genotypes\n";
			return NULL;
		}
		if (os) {
			*os<<"Trio found: Family ID:"<<member->GetFamilyID()<<"\n";
			*os<<"Father:   ";
			father->Report( os );
			*os<<"Mother:   ";
			mother->Report( os );
			*os<<"Child:    ";
			member->Report( os );
		}
		newChild = GetEntry(id);
		newChild->MarkForDeletion( false );
		newChild->SetParents(father, mother);
		newChild->SetStatus( 1 );
		GenotypeData *fakeData 		= newChild->GetGenotypeData();


		uint count = childData->GetGenotypeNumber(0);
	
		int f;
		int m;
		int c1;
		int badData=0;
		for (uint i=0; i<count; i++) {
			f=fatherData->GetGenotypeIndex(i);
			m=motherData->GetGenotypeIndex(i);
			c1=childData->GetGenotypeIndex(i);

			//Let's make sure any parents that are missing data get filtered out as unknown
			if (f * m== 0 || c1==0)
				fakeData->SetGenotype( i, badData);
			//If the parents have the same genotype, we have two options. The children have same, unless both parents are heterozygote
			else if (m==f)
				if (m==2 && c1==1)
					fakeData->SetGenotype(i,3);
				else if (m==2 && c1==3)
					fakeData->SetGenotype(i,1);
				else if (c1 == f)					
					fakeData->SetGenotype( i, c1);
				else {
					if (os)
						*os<<"-\tChild Data("<<i<<") Doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
					fakeData->SetGenotype(i, badData);
				}						
			else {
				//If the parents are seperated by 1, then we have a hetero/homo sitation. Technically, the children will mimic one or the other
				if (abs(m-f) < 2) {
					if (f == c1 || m == c1) {
						if (childData->GetGenotypeIndex(i) == f)
							fakeData->SetGenotype(i, m);
						else
							fakeData->SetGenotype(i, f);
					}	
					else {
						if (os)
							*os<<"-\tChild Data("<<i<<") Doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
						fakeData->SetGenotype(i, badData);
					}						
				}
				//In this case, the parents have different homozygote genotypes. The children must be heterozygous
				else if (c1 == 2) 
					fakeData->SetGenotype(i, 2);
				else {
					if (os)
						*os<<"-\tChild data("<<i<<") doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
					fakeData->SetGenotype(i,badData);
				}
			}
		}
		if (os)	{
			*os<<"New Child:";
			newChild->Report(os);
		}
		//father->MarkForDeletion( true );
		//mother->MarkForDeletion( true );
	} 
	//Indicate that this node should be removed
	else {
		if (os)
			*os<<"Child "<<member->GetFamilyID()<<"x"<<member->GetID()<<" is unaffected. Unable to create a virtual affected sib from an unaffected sibling. All members of the parental group are being marked for deletion. \n";
		member->MarkForDeletion( true );
		father->MarkForDeletion( true );
		mother->MarkForDeletion( true );
		newChild = NULL;
	}
	return newChild;
}

uint FamilyNode::BuildDsps(uint totalIndividualCount) {
	map <string, FamilyMember *>::iterator itr = entries.begin();
	map <string, FamilyMember *>::iterator end = entries.end();

	char pGroupID[1024]; 	
	
	FamilyMember *member;
	Sibship *unit;
	uint count=0;
	uint fatherID=0, motherID=0;

	ExtendedSibship::iterator parentalUnit;

	//Set up the mapping of children to parental groups
	for (; itr != end; itr++) {
		member = itr->second;
		FamilyMember *father = member->GetFather();
		FamilyMember *mother = member->GetMother();
		
		//Make sure there are 2 parents 
		if ((father && mother)) {
			if (father)
				fatherID=atoi(father->GetID());
			if (mother)
				motherID=atoi(mother->GetID());

			sprintf(pGroupID, "%dx%d", fatherID, motherID);
			//string parentalID=pGroupID;
			//check to see if the parental unit new
			parentalUnit = sibships.find(pGroupID);
			if (parentalUnit == sibships.end())	{
				unit = new Sibship(father, mother, totalIndividualCount);
				sibships[pGroupID] = unit;
			}
			else
				unit = parentalUnit->second;
			unit->InsertChild(member);
		}
		else	{
			member->MarkForDeletion( true);
		}
			
	}


	//Iterate over the parental groups and count the number of effective "DSP"s	
	ExtendedSibship::iterator pgend = sibships.end();
	for (ExtendedSibship::iterator i = sibships.begin(); i != pgend; i++) {
		unit=i->second;
		uint affectedCount = unit->affectedCount;
		uint unaffectedCount = unit->unaffectedCount;
		uint localDSPs = affectedCount * unaffectedCount;
		//int localDSPs = unit->affectedCount * unit->unaffectedCount;


		//Let's watch for some trios. Trios will each only be singly represented, so we 
		//don't have to worry about weighting them.
		if (localDSPs == 0) {
			if (unit->affectedCount == 1) {	
				FamilyMember *virtSib = NULL;
				if (expandAllAffecteds || expandTrios)
					virtSib = CreateVirtualSib(unit->children[0]);
				if (virtSib) {
					unit->InsertChild( virtSib);
					count+=2;
					unit->AdjustWeights();
				}	
				else 
					unit->children[0]->MarkForDeletion(true);
			}
			//This is where there is just not a usable set of data
			else {
				for (uint idx=0; idx<unit->children.size(); idx++) 
					unit->children[idx]->MarkForDeletion(true);
			}
		}
		else  {
			//Otherwise, we have a normal parent group. Let's see if we need to duplicate anyone
			//According to Eden Martin, if we bump create virtual sibs for all affected, our power increases 
			//as well. 
			if (expandAllAffecteds) {
				uint childCount = unit->children.size();
				FamilyMember *virtSib = NULL;				//The new child created
				FamilyMember *sib = NULL;					//the current child being considered for expansion
				for (uint i = 0; i<childCount; i++) {
					sib = unit->children[i];
					if (sib->IsAffected()) {
						virtSib = CreateVirtualSib(sib);
						//If the parents don't have data, we can't really do this- but it's not illegal
						if (virtSib) {
							unaffectedCount++;
							unit->InsertChild(virtSib);
						}
					}
				}
			}
			localDSPs = affectedCount * unaffectedCount;
			unit->AdjustWeights();
			count+=localDSPs;
		}
		//Apply the status for each of the units to the families status. This will allow us to scan for counts
	}

	return count;
}


uint FamilyNode::AddGenoData(vector<GenotypeData *> &data,CaseControlStatus pgData[], uint &pgPosition,  bool doScrambleStatus /* = false */) {
	ExtendedSibship::iterator end = sibships.end();
	ExtendedSibship::iterator cur = sibships.begin();

	//We are building our status mask too. So, we need to the starting point for them
	//CaseControlStatus status;			///<Local vector associated with the family
	status.Reset();
	uint pos = data.size();				///<The position within the status bitvector thus far
	//uint childrenConsidered=0;
	
	FamilyMember *children[64];			///<Assuming there are only 64 children to any couple
	uint aIdx=0;						///<Affected index
	uint uIdx=0;						///<Unaffected Index
	uint totalContributed=0;			///<Book keeping so we can tell the caller how many we added
	uint i,n,m;

	dspCount = 0;

	while (cur != end) {
		CaseControlStatus pgStatus;
		Sibship *unit = cur->second;
		cur++;
		//Build up the list of affected children
		uint childCount = unit->children.size();
		assert(childCount<64);

		aIdx=0;
		uIdx=unit->affectedCount;
		for (i=0; i<childCount; i++) {
			
			FamilyMember *child=unit->children[i];
			//We are skipping Deleted and unknown statuses
			if (!child->MarkForDeletion( ) && !child->IsUnknownStatus()) {
				//Sort the children according to status. 
				if (child->IsAffected()) 
					children[aIdx++]=child;
				else
					children[uIdx++]=child;
			}
		}	
		
		//Now, we have a sorted array. We want to pair them up and possibly 
		//shuffle their status now, but only if it makes sense to continue
		if (aIdx > 0 && uIdx != unit->affectedCount)
		{

			//Scramble the status by jumbling the array that is otherwise
			//Sorted according to status
			if (doScrambleStatus) {
				random_shuffle(children, children+uIdx, Utility::Random::globalGenerator);
			}

			//Basic book keeping
			totalContributed+=aIdx*(uIdx-aIdx)*2;
	
			//Let's pair up and add the genetopyes for the ones we are calling "affected"
			for (n=0; n<aIdx; n++) {
				FamilyMember *curChild = children[n];
				for (m=aIdx; m<uIdx; m++) {
					FamilyMember *other = children[m];	
					GenotypeData *gt=curChild->GetGenotypeData(other->GetID());
					assert(gt != NULL);
					//Set the child's new home(s) in the snp repository
					curChild->SetIndividualIdx(pos);	
					curChild->SetEffectiveStatus(FamilyMember::_affectedValue);
					//Set the position in the status mask
					if (pgData)
						pgStatus.SetStatus(pos, true);
					status.SetStatus(pos++, true);
					dspCount++;
					//push the status onto the vector
					data.push_back(gt);
				}
			}
			//Do the same for the "unaffected"
			for (n=unit->affectedCount; n<uIdx; n++) {
				FamilyMember *curChild = children[n];
				for (m=0; m<aIdx; m++) {
					FamilyMember *other =children[m];	
					GenotypeData *gt=curChild->GetGenotypeData(other->GetID());
					assert(gt != NULL);
					//Set the child's new home(s) in the snp repository
					curChild->SetIndividualIdx(pos);
					curChild->SetEffectiveStatus(FamilyMember::_unaffectedValue);

					//Set the position in the status mask
					if (pgData)
						pgStatus.SetStatus(pos, false);
					status.SetStatus(pos++, false);
					//push the status onto the vector
					data.push_back(gt);
				}
			}
		}
		if (pgData)
			pgData[pgPosition++] = pgStatus;
	}
	return totalContributed;
}	

	
}



