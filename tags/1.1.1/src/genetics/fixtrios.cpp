//
// C++ Implementation: fixtrios
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fixtrios.h"

namespace Genetics {

namespace Evaluation {

/**
@brief Removes trio parents from the family and creates the virtual control.
Basically, if the family member is a parent with 1 child, they will be removed from the family.
If the member is a child of a trio group, create the virtual control and add it to the family.
New child will have an id of 25 higher than the original affected child.
*/ 
bool FixTrios::operator()(FamilyMember *member, uint position) {
	familyCount++;
	FamilyMember *father = member->GetFather();
	FamilyMember *mother = member->GetMother();
	
	//If both parents aren't present, then we can't do anything with them
	if ((father == NULL && mother == NULL) || member->MarkForDeletion() )  {
		return false;
	}

	//Verify that the child is the only child by both parents
	if (father->ChildCount() == 1 && mother->ChildCount() == 1) {
		if (father->IsAffected() || mother->IsAffected())	{
			cout<<"Trio found, but one or more parents were found to be affected. Not a valid trio.\n";
			member->MarkForDeletion(true);
			father->MarkForDeletion( true );
			mother->MarkForDeletion( true );
			invalidTrios++;
			return true;
		}
		
		cout<<"Triad found: Family ID:"<<member->GetFamilyID()<<"\n";
		cout<<"Father:   ";
		father->Report( &cout );
		cout<<"Mother:   ";
		mother->Report( &cout );
		cout<<"Child:    ";
		member->Report( &cout );
		if (member->IsAffected()) {
			FamilyMember *newChild;
			//create the copy 
			char id[16];
			sprintf(id, "%d",  atoi(member->GetID())+15);
			
			newChild = family->GetEntry(id);
			newChild->MarkForDeletion( false );
			newChild->SetParents(father, mother);
			newChild->SetStatus( UNAFFECTED );
			GenotypeData *fakeData 		= newChild->GetGenotypeData();
			GenotypeData *childData 	= member->GetGenotypeData();
			GenotypeData *fatherData 	= father->GetGenotypeData();
			GenotypeData *motherData	= mother->GetGenotypeData();

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
						cout<<"-\tChild Data("<<i<<") Doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
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
							cout<<"-\tChild Data("<<i<<") Doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
							fakeData->SetGenotype(i, badData);
						}						
					}
					//In this case, the parents have different homozygote genotypes. The children must be heterozygous
					else if (c1 == 2) 
						fakeData->SetGenotype(i, 2);
					else {
						cout<<"-\tChild data("<<i<<") doesn't match parents: "<<member->GetGtValue(c1)<<" can't be produced from "<<member->GetGtValue(m)<<", "<<member->GetGtValue(f)<<"\n";
						fakeData->SetGenotype(i,badData);
					}
				}
			}
			cout<<"New Child:";
			newChild->Report(&cout);
			father->MarkForDeletion( true );
			mother->MarkForDeletion( true );
			validTrios++;
			childrenCreated++;
			parentsRemoved+=2;
		} 
		//Indicate that this node should be removed
		else {
			cout<<"Child's status is false, family is being marked for removal!\n";
			member->MarkForDeletion(true);
			father->MarkForDeletion( true );
			mother->MarkForDeletion( true );
			invalidTrios++;
			return true;
		}
	}
	return false;
}
			
bool ReportFamilies::operator()(FamilyMember *member, uint position) {
	member->Report(os);
	return true;
}

bool WriteSelectedLoci::VerifyLocus( uint locus ) {
	bool isValid = false;

	if (locus < validLoci.size())
		isValid = validLoci[locus];

	return isValid;
}

bool WriteSelectedLoci::SetupModel( SnpAligned *model) {
	uint totalLoci = model->CountGenotypes();
	uint modelSize = model->GetLabelCount();
	//Set all of the bits to false
	validLoci.reset(totalLoci);
	
	for (uint i=0; i<modelSize; i++) {
		validLoci[model->GetLabel( i )] = true;
	}
	return true;
}

bool WriteFamilies::Evaluate(FamilyMember *member, uint position) {
	bool wasWritten = false;
	string gtValues[] = {"0\t0", "1\t1", "1\t2", "2\t2"};
	GenotypeData *data=member->GetGenotypeData();
	
	uint count =data->GetGenotypeNumber(0);

	if (count > 0) {

	//if (!member->MarkForDeletion()) {
		familyCount++;
		string p1="0";
		string p2="0";
		
		FamilyMember *p=member->GetFather();
		if (p)
			p1=p->GetID();
		p=member->GetMother();
		if (p)
			p2=p->GetID();
		*os<<member->GetFamilyID()<<"\t"<<member->GetID()<<"\t"<<p1<<"\t"<<p2<<"\t0\t0\t0\t0\t0\t";
		if (member->IsAffected())
			*os<<member->_affectedValue<<"\t";
		else
			*os<<member->_unaffectedValue<<"\t";

		for (uint i=0; i<count; i++)  {
			if (VerifyLocus(i)) {
				//cout<<"Writing locus "<<i<<" ("<<gtValues[data->GetGenotypeIndex(i)]<<", "<<data->GetGenotypeIndex(i)<<")\n";
				*os<<gtValues[data->GetGenotypeIndex(i)];
				if (i<count-1)
					*os<<"\t";
			}
		}
		*os<<"\n";
		wasWritten=true;
	}
	return wasWritten;
}



}

}
