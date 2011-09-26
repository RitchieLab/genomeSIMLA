//
// C++ Implementation: pedigreecleaner
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreecleaner.h"
#include <sstream>
#include <iomanip>

namespace MdrPDT {

namespace Validator {

using namespace std;
using namespace Utility;

PedigreeCleaner::PedigreeCleaner(int level, int threshold) :
		level(level), threshold(threshold)	{
}


PedigreeCleaner::~PedigreeCleaner()	{
}

void PedigreeCleaner::Evaluate(PedigreeRepository *repo, std::ostream& os) {
	PedigreeRepository::Iterator itr = repo->GetIterator();
	Pedigree *ped = itr.GetNext();
	int totalErrors = 0;
	errorCounts.clear();
	os<<"----------- Mendelian Error Report ---------------------\n";
	stringstream summary;
	while (ped) {
		totalErrors+=Evaluate(ped, os, summary);
		ped = itr.GetNext();
	}
	os<<summary.str();
	os<<setw(45)<<"Total Genotype Errors: "<<totalErrors<<"\n";
	cout<<setw(45)<<"Total Genotype Errors: "<<totalErrors<<"\n";
}

void PedigreeCleaner::StripLocus(Pedigree *ped, int locus) {
	Pedigree::Iterator itr = ped->GetIterator();
	Individual *ind = itr.GetNext();

	while (ind) {
		ind->SetGenotype(locus, 0);
		ind=itr.GetNext();
	}
}

int PedigreeCleaner::Evaluate(Pedigree* ped, std::ostream& os, std::ostream& summary) {
	Pedigree::Iterator itr = ped->GetIterator();
	Individual *ind = itr.GetNext();
	
	int totalGenotypes = 0;
	int totalErrors = 0;
	
	BitSetType problemLoci;
	if (ind) {
		totalGenotypes = ind->CountGenotypes();
		problemLoci.resize(totalGenotypes, false);
		if (errorCounts.size() == 0)
			errorCounts.resize(totalGenotypes, 0);
		stringstream ss;
		while (ind) {
			totalErrors+=Evaluate(ind, ss, problemLoci);
			ind = itr.GetNext();
		}
		int problemCount = problemLoci.count();
	
		if (problemCount > 0) {
			os<<"Pedigree: "<<ped->ID()<<"\t"<<problemCount<<"/"<<totalGenotypes<<" Affected Loci [ ";
		
			for (int i=0; i<totalGenotypes; i++) {
				if (problemLoci[i]) {
					errorCounts[i]+=1;
					os<<i<<" ";
					if (level>1)
						StripLocus(ped, i);
				}
			}
			os<<"] \n"<<ss.str();
		}
		else 
			os<<"Pedigree: "<<ped->ID()<<"\tNo Errors\n";
		itr.Reset();
		ind=itr.GetNext();
		if (level ==3 && totalErrors > threshold && threshold > 0) {
			os<<"  * Purging "<<ped->GetMemberCount()<<" members due to questionable genotyping information.\n";
		
			while (ind) {	
				ind->DropFromAnalysis(true);
				ind = itr.GetNext();
			}		
			ped->DropFromAnalysis(true);
		}
		else if (level > 1 && problemCount > 0) {
			os<<"  * Stripping "<<problemCount<<" loci from entire pedigree.\n";
		}
	}
	return totalErrors;
}



int PedigreeCleaner::Evaluate(Individual *ind, std::ostream& os, BitSetType& problemLoci) {
	int badData=0;
	if (ind) {
		Individual *mother = ind->GetMother();
		Individual *father = ind->GetFather();

		int missingGts=0;					///<Keep up with the number of empty columns in the child 
		int missingParentData=0;			///<Keep up with the parent's data

		stringstream ss;
		//We can't really tell much about genotypes, if we don't have parents
		if (mother && father) {
			int genotypeCount = ind->CountGenotypes();
			int m=0, f=0, t=0;
			for (int i=0; i<genotypeCount; i++) {
				m=mother->GetGenotype(i);
				f=father->GetGenotype(i);
				t=ind->GetGenotype(i);
				if (t > 0) { 
					//We can't do much if anyone is missing data
					if (f * m > 0) {
						bool valid=true;
						//If the parents have the same genotype, we have two options. 
						//The children have same, unless both parents are heterozygote
						if (m==f)  {
							//If parents are both heterozygous, then we are fine
							if (m==2) 
								valid=true;
							else 
								valid=t==m;
						}
						//In this case, parents are different
						else {
							//If the parents are seperated by 1, then we have a hetero/homo sitation. Technically, the children will mimic one or the other
							if (abs(m-f) == 1) 
								valid=(f == t || m == t);
							//Otherwise, we have two different homozygous parents...can only create heterozyous children
							else 
								valid=t==2;
						}
						//If there is concern for a given locus, we mark that locus as questionable for the pedigree
						if (!valid) {
							problemLoci[i]=true;
							ss<<"("<<i+1<<") "<<t<<"  "<<f<<" x "<<m<<"\t";
							badData++;
						}
					}
					else {
						missingParentData++;
					}
				}
				else {
					missingGts++;
					if (f*m == 0) {
						missingParentData++;
					}
				}
			}
		}

		os<<" - "<<ind->GetID()<<"\t Missing: P)"<<missingParentData<<" C)"<<missingGts<<" ";
		if (badData == 0)
			os<<"\tno mendelian errors found\n";
		else
			os<<"\tErrors: "<<badData<<" \t "<<ss.str()<<"\n";
		
	}
	return badData;
}

}

}
