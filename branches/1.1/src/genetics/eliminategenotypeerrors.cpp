//
// C++ Implementation: fixeliminategenotypeerrors
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "eliminategenotypeerrors.h"

namespace Genetics {

namespace Evaluation {

void FixGenotypeErrors::SetLociCount(uint number) {
	affectedLoci.resize(number, true);

	if (number > affectedCounts.size())
		affectedCounts.resize(number - affectedCounts.size(), 0);
}

void FixGenotypeErrors::Reset() {
	SetLociCount(0);

	//Reset families worst number of genotype errors 
	maxErrorCount = 0;

	totalGenotypes = 0;

	totalPedErrors = 0;

}

void FixGenotypeErrors::Report(ostream *os) {
	*os<<setw(45)<<"Total Mendelian Errors: "<<overallErrors<<endl;
	*os<<setw(45)<<"Level: "<<level;
	switch (level) {
		case 1: 
			cout<<" - No action to be taken\n";
			break;
		case 2:
			cout<<" - Errored Genotypes will be removed\n";
			break;
		case 3:
			cout<<" - Errors will be removed\n";
			cout<<setw(45)<<"Pedigree Threshold: "<<threshold<<" - Any pedigrees that exceed this threshold will be dropped\n";
			break;
		default:
			cout<<" - Unknown level. All errors will be ignored (but reported).\n";
	}
}

bool FixGenotypeErrors::Evaluate(FamilyNode *family, uint position) {
	this->family = family;
	repoCount++;

	Reset();
	
	family->PerformEvaluation(this);
	overallErrors+=totalPedErrors;

	uint count = affectedLoci.size() - affectedLoci.count();

	if (report) 
		*report<<"Family: "<<family->GetID()<<"\t"<<count<<"/"<<totalGenotypes<<" Affected Loci [ ";
	
	for (uint i=0; i< affectedLoci.size(); i++) {
		if (!affectedLoci[i]) {
			affectedCounts[i]=affectedCounts[i]+1;
			if (report)
				*report<<i<<" ";
		}
	}
	if (report)
		*report<<"] ";

	if (level ==3 && totalPedErrors > threshold) {
		ClearBadFamilies clearFam;
		if (report) 
			*report<<" * Purging "<<family->GetMemberCount()<<" members due to questionable genotyping information!";
		family->PerformEvaluation(&clearFam);
		family->MarkForDeletion( true );
	}
	

	//At this point, work through the family, setting all of the members genotypes to unknown where they had to be set
	else if (level > 1 && count > 0) 	{
		StripSibGT stripSibs(&affectedLoci);
		if (report) 
			*report<<"* Stripping "<<count<<" loci!";
		family->PerformEvaluation(&stripSibs);
	}

	if (report) 
		*report<<"\n";

	return false;
}

string FixGenotypeErrors::HandleLoci(uint loci, FamilyMember *member, uint father, uint mother, uint child) {
	//int errorIdx=0;
	char values[128];
	sprintf(values, "(%d) %s   %s x %s\t", loci, member->GetGtValue(child).c_str(), member->GetGtValue(father).c_str(), member->GetGtValue(mother).c_str());
	
	//if (level > 1)
	//	member->GetGenotypeData()->ZeroGenotype( loci);

	badData++;
	affectedLoci[loci]=false;
	//affectedCounts[loci]=affectedCounts[loci]+1;
	return values;
}



/**
@brief Removes trio parents from the family and creates the virtual control.
Basically, if the family member is a parent with 1 child, they will be removed from the family.
If the member is a child of a trio group, create the virtual control and add it to the family.
New child will have an id of 25 higher than the original affected child.
*/ 
bool FixGenotypeErrors::Evaluate(FamilyMember *member, uint position) {
	badData=0;
	familyCount++;
	FamilyMember *father = member->GetFather();
	FamilyMember *mother = member->GetMother();

	GenotypeData *childData 	= member->GetGenotypeData();
	uint count = childData->GetGenotypeNumber(0);
	SetLociCount(count);
	if (totalGenotypes < count)
		totalGenotypes = count;
	
	//If both parents aren't present, then we can't do anything with them. 
	//This basically means two things: parents or discordant sibships. 
	//If they don't have children, we treat them as discordant sibships. If this 
	//is the case, we want to make sure they aren't set to be ignored
	if (father == NULL && mother == NULL)  {
		if (member->ChildCount() == 0)	 {
			*report<<" - "<<member->GetID()<<"\t Missing: P) (no data) C) (all unchanged) \n";
			member->MarkForDeletion(false);
		}
		return false;
	}

	//*report<<"Father:   ";
	//father->Report( report );
	//*report<<"Mother:   ";
	//mother->Report( report );
	//*report<<"Child:    ";
	//member->Report( report );


	GenotypeData *fatherData = NULL;
	if (father)
		fatherData = father->GetGenotypeData();
		
	GenotypeData *motherData = NULL;
	if (mother)
		motherData = mother->GetGenotypeData();

	int f, m, c1;
	int missingGts=0;					///<Keep up with the number of empty columns in the child 
	int missingParentData=0;
	
	
	string indReport="";

	for (uint i=0; i<count; i++) {
		f=0;
		if (fatherData)
			f=fatherData->GetGenotypeIndex(i);

		m=0;
		if (motherData)
			m=motherData->GetGenotypeIndex(i);

		c1=childData->GetGenotypeIndex(i);
			
		if (c1 == 0) {
			//affectedLoci[i]=false;
			missingGts++;
			if (m*f == 0){
				//This case is acceptable- since we can't make any decisions about the kid's heritage
				missingParentData++;
			}
		}
		else if (m+f == 0){
			//This case is acceptable- since we can't make any decisions about the kid's heritage
			missingParentData++;
		}
		//If one parent is missing data, we need to check for Heterozygous conditions
		else if (f * m== 0){
			missingParentData++;
			if (f > 0) {
				//If the difference of the two is more than one, then we have impossible homozygous conditions
				if (abs(f-c1)> 1) {										indReport+=HandleLoci(i, member, f, m, c1);  }
			}	
			else {
				if (abs(m-c1) > 1) {									indReport+=HandleLoci(i, member, f, m, c1);  }
			}
		}
		//If the parents have the same genotype, we have two options. The children have same, unless both parents are heterozygote
		else if (m==f) {
			if (! ((m==2 && c1==1)|| (m==2 && c1==3) || (c1 == f))	) {	indReport+=HandleLoci(i, member, f, m, c1);  }
		}
		//If the parents are seperated by 1, then we have a hetero/homo sitation. Technically, the children will mimic one or the other
		else if ((abs(m-f) < 2)) {
			if (!(f == c1 || m == c1)) {								indReport+=HandleLoci(i, member, f, m, c1);  }
		}
		//In this case, the parents have different homozygote genotypes. The children must be heterozygous
		else if (c1 != 2)  {									 		indReport+=HandleLoci(i, member, f, m, c1);  }
	}

	if (report) {
		*report<<" - "<<member->GetID()<<"\t Missing: P)"<<missingParentData<<" C)"<<missingGts<<" ";
		if (indReport.size() == 0)
			*report<<"\tno mendelian errors found\n";
		else
			*report<<"\tErrors: "<<badData<<" \t "<<indReport<<"\n";
	}	
	totalPedErrors+=badData;

	return false;
}

void StripBadLoci::Init() {
	if (report.size() < affCounts->size()) 
		for (uint i=report.size(); i<affCounts->size(); i++)
			report.push_back(0);
}

/**
@brief Removes trio parents from the family and creates the virtual control.
Basically, if the family member is a parent with 1 child, they will be removed from the family.
If the member is a child of a trio group, create the virtual control and add it to the family.
New child will have an id of 25 higher than the original affected child.
*/ 
bool StripSibGT::Evaluate(FamilyMember *member, uint position) {
	familyCount++;
	FamilyMember *father = member->GetFather();
	FamilyMember *mother = member->GetMother();

	//Skip Parents
	if (father == NULL || mother == NULL)  {
		return false;
	}

	uint count = member->GetGenotypeData()->GetGenotypeNumber(0);
		
	for (uint i=0; i<count; i++) {
		if (! affectedLoci->test(i)) {
			member->GetGenotypeData()->ZeroGenotype(i);
		}
	}
	
	return false;
}



/**
 * @brief Identifies any invalid genotypes for each child (with parent information) and zeros them out
 */
bool StripBadLoci::Evaluate(FamilyMember *individual, uint position) {
	Init();
	uint count = affCounts->size();
	GenotypeData *childData 	= individual->GetGenotypeData();

	bool success= false;	

	for (uint i=0; i<count; i++) {
		if (affCounts->at(i) > threshold && childData->GetGenotypeNumber(i) != 0) {
			report[i] = report[i]++;
			if (mode == 3) {
				childData->ZeroGenotype(i);
				//childData->SetGenotype(i, (uint)0);
//				cout<<"Setting Locus #"<<i<<" to all zeros\n";
			}
			success=true;
		}
	}
	return success;
}
	 
/**
 * Generate complete sumary
 */
void StripBadLoci::Report(ostream *os) {
	uint count = report.size();

	*os<<"\n\n----------------------------Loci Cleaning\n";
	*os<<"Stripping loci exceeding threshold of "<<threshold<<"\n";
	for (uint i=0; i<count; i++) {
		if (report[i]>0)
			if (mode == 3)
				*os<<i<<": "<<report[i]<<" individuals cleared\n";
			else
				*os<<i<<": "<<report[i]<<" would have been cleared\n";
	}
}

/**
 * @brief used to reset the affected loci vector
 */
void StripBadLoci::Reset() {
	uint count = report.size();

	for (uint i=0; i<count; i++) {
		report[i]=0;
	}
}


}

}
