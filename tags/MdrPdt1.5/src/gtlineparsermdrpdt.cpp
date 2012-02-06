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
#include "sibshipfoldproduction.h"
#include <iomanip>

namespace MDR {
bool GtLineParserMdrPdt::WriteFoldDetails = false;

void GtLineParserMdrPdt::ExcludePedigrees(vector<string> pedigreeList) {
	pedExclusionList = pedigreeList;
	
}

void GtLineParserMdrPdt::ExcludeLoci(vector<string> lociExclusions) {
	lociExclusionList = lociExclusions;
}

void GtLineParserMdrPdt::PostParse() {
	uint count = pedExclusionList.size();
	string pedigreeID;
	//Let's drop any families that were added to the list of exclusions
	if (count > 0) {
		for (uint i=0; i< count; i++) {
			FamilyNode *fam = famRepo->GetNode(pedExclusionList[i].c_str());
			if (fam) {
				fam->MarkForDeletion(true);
				*pedStream<<"Removing family: "<<fam->GetID()<<". "<<fam->GetMemberCount()<<" individuals ignored\n";
			}
		}
	}

	count = lociExclusionList.size();
	snpExclusions.resize(lociObserved, false);
	snpExclusions.reset();

	for (uint i=0 ; i<count; i++) {
		int locus = atoi(lociExclusionList[i].c_str());
		if (locus > 0)
			snpExclusions[locus - 1] = true;
	}



	SetupChildren(false);

}



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

	if (lociObserved == 0)
		lociObserved = idx;
	else if (lociObserved != idx) {
		cout<<"Line #"<<val<<" wasn't the same size as previous lines. Please check that there is nothing wrong with the dataset\n";
		abort();
	}
		

	//We don't want the ones that are dropped in from reference to get written so we mark them when we have data for them
	if (!child->IsUnknownStatus())
		child->MarkForDeletion( false );

	return true;
}

void GtLineParserMdrPdt::ReportOnIndividualInUse() {
	if (pedStream) {
		FamilyRepository::Iterator i = famRepo->GetIterator();
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
	uint ignoredSnps = snpExclusions.count();

	uint width=45;
	//CaseControlStatus st=statusMask->CombinedStatus();
	os<<setw(width)<<right<<"Input File Type: "<<"Pedigree Format\n";
	if (FamilyNode::expandAllAffecteds) {
		os<<setw(width)<<right<<" "<<"* Where adequate parental data was avaialble, \n";
		os<<setw(width)<<right<<" "<<"  virtual sibs were created to pair with affected siblings.\n";
	}		
	else if (FamilyNode::expandTrios) {
		os<<setw(width)<<right<<" "<<"* Trios were expanded to valid\n";
		os<<setw(width)<<right<<" "<<"  DSPs where parental data was present\n";
	}
	if (pedExclusionList.size() > 0) {
		os<<setw(width)<<right<<"Pedigrees to be ignored: [ ";
		copy(pedExclusionList.begin(), pedExclusionList.end(), ostream_iterator<string>(os, " "));
		os<<" ]\n";
	}

	os<<setw(width)<<right<<"Value Denoting Affected: "<<FamilyMember::_affectedValue<<"\n";
	os<<setw(width)<<right<<"Value Denoting Unaffected: "<<FamilyMember::_unaffectedValue<<"\n";
	os<<setw(width)<<right<<"Total Individuals Seen: "<<totalIndividualsSeen<<"\n";
	os<<setw(width)<<right<<"Total individuals to be used in analysis: "<<includedChildren.size()<<"\n";

	os<<setw(width+5)<<right<<"Affected Individuals: "<<status.affected.count()<<"\n";
	os<<setw(width+5)<<right<<"Unaffected Individuals: "<<status.unaffected.count()<<"\n";

	if (ignoredSnps > 0) {	
		os<<setw(width)<<right<<"Number of SNPS encountered: "<<GetSnpCount()<<"\n";
		os<<setw(width)<<right<<"Snps Ignored: "<<ignoredSnps<<"[ ";
		copy(lociExclusionList.begin(), lociExclusionList.end(), ostream_iterator<string>(os, " "));
		os<<"] \n";
	}
	os<<setw(width)<<right<<"Number of SNPS to be analyzed: "<<(GetSnpCount() - ignoredSnps)<<"\n";

	if (pedStream) {
		famRepo->GenerateReport(pedStream);
		*pedStream<<"\t* Indicates that the child at that ID is affected\n";
		*pedStream<<"\t? Indicates that the child at that ID has unknown status\n\n";
	}
}

//We are going to move the role of status assignment over to the fold generation
void GtLineParserMdrPdt::SetupChildren(bool doRandomize) {
	CaseControlStatus *pgStat=NULL;

	FamilyRepository::Iterator i = famRepo->GetIterator();
	FamilyNode *family = i.GetNext();
	includedChildren.clear();

	uint totalInUse=0;
	if (!doRandomize) {
		while (family) {
			if (!family->MarkForDeletion())
				if (!family->BuildDsps(totalIndividualsSeen) > 0)
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
			totalInUse+=family->AddGenoData(includedChildren, pgStat, pgStatArrayCount, doRandomize);
		}
		family=i.GetNext();
	}	

	//We don't know how many children we have until we finish the AddGenoData loop. 
	//--This isn't entirely true. After we run this once, this value won't change. This
	//could be a source for improvement in the future
	uint childCount = includedChildren.size();

	if (childCount > 0) {


		//Just for testing purposes, let's create a Sibship Fold Production object and report on it's findings
		SibshipFoldProduction foldProduction(foldCount,famRepo, includedChildren.size()); 
		if (folds)
			delete[] folds;
		folds = foldProduction.BuildXVFolds(childCount);
#ifdef CROSS_VALIDATION
	if (WriteFoldDetails && pedStream) {
		foldProduction.GenerateReport(pedStream);
	}
#endif
	
		//Set up the overall status
		for (uint i=0; i<foldCount;i++)  {
			CaseControlStatus status = folds[i].GetOverallStatus();
			overAllStatus.AppendStatus(status);
		}
	
		if (pgStat) {
			for (uint i=0; i<pgStatArrayCount; i++) {
				pgStat[i].Resize(childCount);
			}
		}
	}
	else 
		cout<<"No valid entries were found in the data. \n";
} 



SnpAligned *GtLineParserMdrPdt::GetSnp(uint snpID) {
	SnpAligned *snp=NULL;
	//This function assumes that all children that are to be included are stored in some structure somehow
	uint childCount=includedChildren.size();

	//If we have nothing inside our buffer, then we can't do much else!
	if (childCount == 0)
		return snp;

	if (snpID < includedChildren[0]->GetGenotypeNumber(0)) {
		//Acquire an empty snp from the pool
		snp=pool->GetSnp(genoCount, childCount);
		snp->Reset();
		//If this is to be excluded, let's let it end up with no data
		if (!snpExclusions[snpID]) {
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
		}
	}
	
	return snp;		
}


}
