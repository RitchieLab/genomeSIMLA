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
#include "genetics/fixtrios.h"
#include "strings.h"

namespace MDR {
bool GtLineParserMdrPdt::WriteFoldDetails = false;

void GtLineParserMdrPdt::ExcludePedigrees(vector<string> pedigreeList) {
	pedExclusionList = pedigreeList;
	
}
void GtLineParserMdrPdt::WriteSaScript(const char *filename, vector<SnpAligned*> loci, 
		const char *targetPed, bool withInteraction) {
	ofstream os(filename, ios::out);

	if (!os.is_open()){
		cout<<"Unable to open file: "<<filename<<" to write sas script\n";
		abort();
	}

	uint modelSize = loci.size();
	vector<string> varComponents;				///<These are used to construct the interaction variables

	os<<"PROC IMPORT OUT = WORK.PEDTEST\n";
	os<<"            DATAFILE=\""<<targetPed<<"\"\n";
	os<<"            DBMS=TAB REPLACE;\n";
	os<<"     GETNAMES=YES;\n";
	os<<"     DATAROW=2;\n";
	os<<"     GUESSINGROWS=20;\n";
	os<<"RUN;\n\n";
	os<<"data work.ped; set work.pedtest;\n";
	os<<"sibship=compress(famid||fath||moth);\n";
	

	//Set up the missing data 
	stringstream columns;
	stringstream dbVars;
	stringstream products;
	stringstream hrResults;

	stringstream interactionTerm;
	

	for (uint i = 0; i<modelSize; i++)  {
		os<<"\n";
		SnpAligned *model = loci[i];
		BitSetType highRiskCells = model->GetHrCells();
//cout<<"!!-->"<<model->GetLabel()<<" : "<<highRiskCells<<"\n";
		uint genotypeCount = model->CountGenotypes();
		string prd, prefix;
		char buff[1024];

		sprintf(buff, "M%d_", i+1);
		prefix = buff;
		varComponents.push_back(buff);

		sprintf(buff, "%sx", prefix.c_str());
		prd = buff;

		interactionTerm<<prefix;

		os<<"if "<<prefix<<"1=0 and "<<prefix<<"2=0 then "<<prefix<<"x=.;\n";


		if (i>0) {
			hrResults<<" ";
			dbVars<<" ";
		}
		dbVars<<"db"<<prd<<" ";
		hrResults<<prd;
		
		for (uint gt=1; gt<genotypeCount; gt++) {
			//Convert "1/1" into "M1_1" M1_2=1
			string curGenotypes = model->GetGenotypeLabel(gt);
			string v1, v2;

			sprintf(buff, "%s1", prefix.c_str());
			v1 = buff;

			sprintf(buff, "%s2", prefix.c_str());
			v2 = buff;

			uint isHighRisk = highRiskCells[gt];

	
			sprintf(buff, "if %s=%c and %s=%c then %s=%d;", v1.c_str(), curGenotypes.at(0), v2.c_str(), curGenotypes.at(2), prd.c_str(), isHighRisk);
			os<<buff<<"\n";
			
			columns<<v1<<" "<<v2<<" ";
		}	
	}


	//We need to keep up with the size of the models we are describing for order
	stringstream secondaryVariables;
	stringstream dbSecondaryVariables;
	stringstream basicVariables;
	stringstream dbVariables;
	stringstream compoundVariableL;
	stringstream compoundVariableR;
	
	os<<"\n";

	//If we want this to be able to support >3 loci, we need to make it a recursive function
	//Construct the interactions
	for (size_t i=0; i<modelSize; i++) {
		string var=varComponents[i];
		basicVariables<<var<<"x ";
		dbVariables<<"db"<<var<<"x ";

		if (withInteraction) {
			compoundVariableL<<var;
			if (i>0)
				compoundVariableR<<"*";
			compoundVariableR<<var<<"x";
			
			for (size_t j=i+1; j<modelSize; j++) {
				string var2=varComponents[j];
				os<<var<<var2<<"x="<<var<<"x*"<<var2<<"x;\n";			
				secondaryVariables<<var<<var2<<"x ";
				dbSecondaryVariables<<"db"<<var<<var2<<"x ";
			}
		}
	}
	
	if (withInteraction) {
		size_t depth =  modelSize-1;
		basicVariables<<secondaryVariables.str()<<" ";
		dbVariables<<dbSecondaryVariables.str()<<" ";

		if (withInteraction && modelSize > 2) {
			basicVariables<<compoundVariableL.str()<<"x ";	
			dbVariables<<" db"<<compoundVariableL.str()<<"x ";
			
			os<<compoundVariableL.str()<<"x="<<compoundVariableR.str()<<";\n";
		}
	}
/*
	if (withInteraction) {
		hrResults
		interactionTerm<<"x";
		hrResults<<" "<<interactionTerm.str();

		//I still need to work out how to do 3 way interactions;
		os<<interactionTerm.str()<<"="<<products.str()<<";\n";
		
		dbVars<<"db"<<interactionTerm.str();

	}
	*/	

	os<<"\nproc sort data=ped; \n";
	os<<"by sibship;\n";

	os<<"data aff;\n";
	os<<"	set ped;\n";
	os<<"	if aff_status=2;\n";
	os<<"proc sort data=aff;\n";
	os<<"	by sibship; \n";
	os<<"data aff_n (keep=famid sibship affsib_n);\n";
	os<<"	set aff;\n";
	os<<"	by sibship;\n";
	os<<"	retain affsib_n;\n";
	os<<"	if first.sibship then affsib_n=0;\n";
	os<<"	affsib_n=affsib_n+1;\n";
	os<<"	if last.sibship then output;\n";

	os<<"data unaff;\n";
	os<<"	set ped;\n";
	os<<"	if aff_status=1;\n";
	os<<"proc sort data=unaff;\n";
	os<<"	by sibship; \n";
	os<<"data unaff_n (keep=famid sibship unaffsib_n);\n";
	os<<"	set unaff;\n";
	os<<"	by sibship;\n";
	os<<"	retain unaffsib_n;\n";
	os<<"	if first.sibship then unaffsib_n=0;\n";
	os<<"	unaffsib_n=unaffsib_n+1;\n";
	os<<"	if last.sibship then output;\n";

	os<<"data dsp;\n";
	os<<"	merge aff_n (in=a) unaff_n (in=b);\n";
	os<<"	by sibship;\n";
	os<<"	if a and b;\n";
	os<<"	dsp_n=affsib_n+unaffsib_n;\n";

	os<<"proc sort data=dsp;\n";
	os<<"	by sibship;\n";

	os<<"data ped_dsp;\n";
	os<<"	merge ped (in=a) dsp (in=b);\n";
	os<<"	by sibship;\n";
	os<<"	if a and b;\n";
	os<<"	if aff_status=2 then case=1;\n";
	os<<"	else if aff_status=1 then case=0;\n";
	os<<"	else case=case;\n";
	os<<	"time=2-case;\n";

	os<<"proc sort data=ped_dsp;\n";
	os<<"	by sibship;\n";

	os<<"proc phreg data=ped_dsp outest=est4 nosummary;\n";
	os<<"	where M1_x ne . and M2_x ne .;\n";
	os<<"	model time*case(0)= "<<basicVariables.str()<<" /rl ties=discrete; \n";
	os<<"	strata sibship;                                         \n";
	os<<"	output out=resids dfbeta="<<dbVariables.str()<<" /order=data;\n";
	os<<"	id sibship;\n";
	os<<"proc sort data=resids;\n";
	os<<"	by sibship;\n";
	os<<"proc means data=resids noprint;\n";
	os<<"	by sibship;\n";
	os<<"	var "<<dbVariables.str()<<"; \n";
	os<<"	output out=out4 sum="<<dbVariables.str()<<"; \n";
	os<<"proc iml;\n";
	os<<"	use out4;\n";
	os<<"	read all var {"<<dbVariables.str()<<"} into r; \n";
	os<<"	rvar=r` * r;\n";
	os<<"	print rvar [rowname={ "<<dbVariables.str()<<"} colname={"<<dbVariables.str()<<"}]; \n";
	os<<"	create est1 from rvar [rowname={"<<dbVariables.str()<<"} colname={"<<dbVariables.str()<<"}]; \n";
	os<<"	quit;\n";


	os.close();
}

void GtLineParserMdrPdt::WritePedFile(const char *filename) {
	ofstream ped(filename, ios::out);
	WriteFamilies famWriter (&ped, false);
	famRepo->PerformEvaluation(&famWriter);
	ped.close();
	
}
void GtLineParserMdrPdt::WriteModelLociOnly(const char *modelLoci, const char *filename) {
	BitSetType validLoci(lociObserved, false);

	char *start = (char *)modelLoci;					//Start of a given number in string
	char *lastChar = start + strlen(start);
	char *end = strchr(start, 'x');	//End position of the number in the string
	
	uint modelID;
	char snpID[64];
	
	bool doContinue = true;	
	uint modelSize = 0;
	//Build the model 1 step at a time
	while (doContinue) {
		uint len = lastChar - start;
		if (end)
			len=end-start;
		strncpy(snpID, start, len);
		snpID[len]='\0';
		modelID = atoi(snpID);

		//We are expecting 1 based indices (equal to or less than the number of snps)
		assert(lociObserved >= modelID && modelID > 0);
		validLoci[modelID - 1] = true;
		modelSize++;

		doContinue = end != NULL;
		if (doContinue) {
			start = end + 1;
			end = strchr(start, 'x');
		}
	}
	ofstream os(filename, ios::out);

	//Write header required by sas
	os<<"famid\tind\tfath\tmoth\toff\tpat_sib\tmat_sib\tsex\tproband\taff_status"; 	//M1_1	M1_2	M2_1	M2_2
	for (uint i=0; i<modelSize; i++) 
		for (uint j=0; j<2; j++) 
			os<<"\tM"<<i+1<<"_"<<j+1;
	//End of header
	os<<"\n";
	WriteSelectedLoci wr(&os, validLoci, false);

	famRepo->PerformEvaluation(&wr);
	os.close();
}

void GtLineParserMdrPdt::ExcludeLoci(vector<string> lociExclusions) {
	lociExclusionList = lociExclusions;
}

void GtLineParserMdrPdt::PostParse() {


	uint evalCount = postEvalMethods.size();

	for (uint i=0; i<evalCount; i++)
		famRepo->PerformEvaluation(postEvalMethods[i]);


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
	string data=line;

	vector<string> tokens;
	Utility::TokenizeString(line, tokens, " \t,\r");
	vector<string>::iterator i = tokens.begin();
	vector<string>::iterator end = tokens.end();

	
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
	vector<string> tokens;
	Utility::TokenizeString(line, tokens, " \t,\r");
	vector<string>::iterator i = tokens.begin();
	vector<string>::iterator end = tokens.end();


	string convertedLine=line;

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
		//We should skip families that are deleted- if we are doing a random run, we don't want to build DSPs
		if (!family->MarkForDeletion()) {
			totalInUse+=family->AddGenoData(includedChildren, pgStat, pgStatArrayCount, doRandomize);

		}
		family=i.GetNext();
	}	

/*	CaseControlStatus temp;
	for (size_t i=0; i<totalInUse; i++) {
		temp.ForceAppend( pgStat[i]);
	}
	cout<<"Setup Children- PG Status ";
	if (doRandomize) 
		cout<<" (Random) ";
	cout<<"\n";
	cout<<"\tAffected   :"<<temp.affected<<"\n";
	cout<<"\tUnaffected :"<<temp.unaffected<<"\n";
*/

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
	if (!doRandomize && WriteFoldDetails && pedStream) {
		foldProduction.GenerateReport(pedStream);
	}
	
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
		cout<<"No valid entries were found in the data. Do status values match the configuration file? \n\n";

/*	cout<<"Setup Children - Overall Status";
	if (doRandomize) 
		cout<<" (Random) ";
	cout<<"\n";
	cout<<"\tAffected   :"<<overAllStatus.affected<<"\n";
	cout<<"\tUnaffected :"<<overAllStatus.unaffected<<"\n";
	*/
} 


void GtLineParserMdrPdt::AddPostEvaluation(FamilyRepoEvaluation *eval) {
	postEvalMethods.push_back(eval);
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
				size_t position = lostData.find_first();
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
