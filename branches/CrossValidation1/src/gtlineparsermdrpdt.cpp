//
// C++ Implementation: gtlineparsermdrpdt
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gtlineparsermdrpdt.h"
#include "genetics/familymember.h"
#include <iomanip>

namespace ESE {

/*
 * This is called in order to set up labels other than line number (+1) for the snps
 */
bool GtLineParserMdrPdt::InitData(const char *line) {
	boost::char_separator<char> sep(" \t,", "", boost::drop_empty_tokens);;
	string data=line;
	strtokenizer tok(data, sep);						
	strtokenizer::iterator i = tok.begin();
	strtokenizer::iterator end = tok.end();
	
	int snpCount=0;
	
	
	labels.clear();

	for (; i != end; i++, snpCount++) {
		if (readLabel)	
			labels.push_back(SnpAligned::Label(snpCount, atoi((*i).c_str())));
		else {
			labels.push_back(SnpAligned::Label(snpCount, snpCount+1));
		}
	}	
	return readLabel;
}


bool GtLineParserMdrPdt::ParseLine( const char *line, uint val) {
	//If this is the first column, we need to populate our members
	if (val==0)		{
		if (InitData(line)) {
			return false;
		}
	}
	
	if (strlen(line) < 1) {
		return false;
	}
	totalIndividualsSeen++;

	string convertedLine=line;
	boost::char_separator<char> sep(" \t,", "", boost::drop_empty_tokens);
	strtokenizer tok(convertedLine, sep);						
	strtokenizer::iterator i = tok.begin();
	strtokenizer::iterator end = tok.end();


	string familyID=*i++;
	string indID=*i++;
	string dad=*i++;
	string mom=*i++;
	string firstOffspringID = *i++;
	string nextPatSibID = *i++;
	string nextMatSibID = *i++;
	string gender = *i++;
	string probandStatus = *i++;
	string status= *i++;
	

	FamilyNode 	*family=famRepo->GetNode(familyID.c_str());
	FamilyMember *father = family->GetEntry(dad.c_str());
	FamilyMember *mother = family->GetEntry( mom.c_str());
	FamilyMember *child = family->GetEntry( indID.c_str());
	
	//This is not preferable...will have to hold onto twice the data for now:(
	GenotypeData *childData=child->GetGenotypeData();

	child->SetParents(father, mother);
	child->SetStatus(atoi(status.c_str()));

	uint idx=0;
	string genotype;

	string newValue;
		
	for (; i != end; i++) {	
		genotype=*i+" ";
	
		//Currently, we get these for junk whitespace. We need to clean the iterator
		if (++i == end) {
			cout<<"Malformed pedigree file encountered "<<val<<" - "<<idx<<" ("<<genotype<<"): \n";
			break;
		}
		genotype+=*(i);
		//Set up the child's data- we might change the triad repair to use the SnpAligned data...but for now, this is adequate
		childData->SetGenotype(idx++, genotype.c_str());

	}
	//We don't want the ones that are dropped in from reference to get written so we mark them when we have data for them
	if (!child->IsUnknownStatus())
		child->MarkForDeletion( false );

	return true;
}

void GtLineParserMdrPdt::ReportOnIndividualInUse() {
	if (pedStream) {
		FamilyRepository::Iterator i = famRepo->FirstFamily();
		FamilyNode *family = i.GetNext();
	
		//iterate over each family and let them contribute their portion to the array of statuses
		while (family) {
			if (!family->MarkForDeletion()) {
				family->GenerateReport(pedStream);
			}
			family=i.GetNext();
		}
	}	
}

void GtLineParserMdrPdt::GenerateReport(ostream &os) {
	CaseControlStatus status;
	GetStatusMask(&status);
	uint width=45;
	//CaseControlStatus st=statusMask->CombinedStatus();
	os<<setw(width)<<right<<"Input File Type: "<<"Pedigree Format\n";
	if (doCreateVirtualSibs) {
		os<<setw(width)<<right<<" "<<"* Trios were expanded to valid\n";
		os<<setw(width)<<right<<" "<<"  DSPs where parental data was present\n";
	}
	os<<setw(width)<<right<<"Value Denoting Affected: "<<FamilyMember::_affectedValue<<"\n";
	os<<setw(width)<<right<<"Value Denoting Unaffected: "<<FamilyMember::_unaffectedValue<<"\n";
	os<<setw(width)<<right<<"Total Individuals Seen: "<<totalIndividualsSeen<<"\n";
	os<<setw(width)<<right<<"Total individuals to be used in analysis: "<<includedChildren.size()<<"\n";

	os<<setw(width+5)<<right<<"Affected Individuals: "<<status.affected.count()<<"\n";
	os<<setw(width+5)<<right<<"Unaffected Individuals: "<<status.unaffected.count()<<"\n";
	os<<setw(width)<<right<<"Number of SNPS per individual: "<<GetSnpCount()<<"\n";

	if (pedStream) {
		famRepo->GenerateReport(pedStream);
		*pedStream<<"\t* Indicates that the child at that ID is affected\n";
		*pedStream<<"\t? Indicates that the child at that ID has unknown status\n\n";
	}

}


void GtLineParserMdrPdt::SetupChildren(bool doRandomize) {
	CaseControlStatus *pgStat=NULL;

	FamilyRepository::Iterator i = famRepo->FirstFamily();
	FamilyNode *family = i.GetNext();
	includedChildren.clear();
	//statusArray.clear();
	if (statArray == NULL) 
		statArray = new CaseControlStatus[famRepo->GetFamilyCount()];
	statArrayCount = 0;

	uint totalInUse=0;
	if (!doRandomize) {
		while (family) {
			if (!family->MarkForDeletion())
				if (!family->BuildDsps(totalIndividualsSeen, doCreateVirtualSibs) > 0)
					family->MarkForDeletion(true);
			family = i.GetNext();
		}
		if (pgStatArray)
			delete[] pgStatArray;
		pgStatArray = new CaseControlStatus[famRepo->CountSibships()];
		pgStat=pgStatArray;
		pgStatArrayCount = 0;
	}

	i.Reset();
	family = i.GetNext();

	//iterate over each family and let them contribute their portion to the array of statuses
	while (family) {
		//We should skip families that are deleted- if we are doing a radom run, we don't want to build DSPs
		if (!family->MarkForDeletion()) {
			totalInUse+=family->AddGenoData(includedChildren, statArray, pgStat, statArrayCount, pgStatArrayCount, doRandomize);
		}
		family=i.GetNext();
	}	

	//We don't know how many children we have until we finish the AddGenoData loop. 
	//--This isn't entirely true. After we run this once, this value won't change. This
	//could be a source for improvement in the future
	uint childCount = includedChildren.size();

	//Loop over each item in the array and resize so that we have homogenous genotypes
	for (uint i=0; i<statArrayCount; i++) {
		statArray[i].Resize(childCount);	
		overAllStatus.AppendStatus(statArray[i]);
	}

	if (pgStat) {
		for (uint i=0; i<pgStatArrayCount; i++) {
			pgStat[i].Resize(childCount);
		}
	}
} 



SnpAligned *GtLineParserMdrPdt::GetSnp(uint snpID) {
	SnpAligned *snp=NULL;
	//This function assumes that all children that are to be included are stored in some structure somehow
	uint childCount=includedChildren.size();

	if (snpID < includedChildren[0]->GetGenotypeNumber(0)) {
		//Acquire an empty snp from the pool
		snp=pool->GetSnp(genoCount, childCount);
		//Set the values according the genotypes that are described in our little vector
		for (uint i=0; i<childCount; i++) {
			snp->SetValue(i, includedChildren[i]->GetGenotypeIndex(snpID));
		}
		//Set up the textual details necessary for reporting
		snp->SetLabel(labels[snpID]);
		snp->SetupGenotypeLabels();
		BitSetType lostData=snp->Verify();
		if (lostData.count() > 0) {
			cout<<"Unable to verify snp: "<<snp->GetLabel()<<"\nThe following individuals should have had the following GTs\n";
			uint position = lostData.find_first();
			while (position != BitSetType::npos) {
				cout<<"\t"<<position<<" "<<includedChildren[position]->GetGenotypeIndex(snpID)<<"\n";
				position = lostData.find_next(position);
			}
		}
		//assert(snp->Verify());
		//BitSetType &missing = snp->GetGenotype( 0 );
		//missing=~missing; 
	}
	
	return snp;		
}


}
