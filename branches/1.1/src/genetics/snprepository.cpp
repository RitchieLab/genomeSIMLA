//
// C++ Implementation: snprepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snprepository.h"
#include "utility/utility.h"
#include <fstream>
#include <iostream>
#include "binarygenoparser.h"
#include "snpverificationmethod.h"


namespace Genetics {

using namespace std;

SnpRepository::~SnpRepository()
{
	PurgeSnps();
	pool->Release();
}

void SnpRepository::Resize(int count) {
	int curSize=snps.size();
	SnpAligned *p = NULL;

	if (curSize <= count - 1) {
		snps.resize(count);

		for (int i=count-1; i>curSize-1; i--) 
			snps[i]=NULL;
	}
	else	{
		for (int i=curSize; i<count; i++)	{
			p = snps[i+1];
			pool->ReleaseSnp(p);
		}
		snps.resize(count);
	}
}



int SnpRepository::CountLines(const char *filename)	{
	int linecount = 0;
	ifstream fin(filename, ios_base::in);
	char line[MAX_LINE_LENGTH];
	if (fin.eof() || fin.fail() || fin.bad() ) {
		cout << "Unable to open file, " << filename << ". \n";
		return false;
	}
	
	while (!fin.eof() ){
		fin.getline(line, MAX_LINE_LENGTH);
		if (strlen(line) > 0)
			linecount++;
	}
	fin.close();
	return linecount;
}

void SnpRepository::PurgeSnps() {
	SnpAligned *p=NULL;
	//cout<<"Purging SNPS( "<<snps.size()<<", "<<curSnp<<", "<<snpCount<<")\n";
	if (snps.size() > 0) {

		for (uint i=0; i<curSnp; i++) {
			//cout<<"Dropping SNP("<<i<<"): ";
			p=snps[i];
			if (p)	{
				//cout<<p->GetID()<<"\n";
				pool->ReleaseSnp(p);
			} else { 
				//cout<<" NULL \n";
			}
		}
	}
	snpLookup.clear();

	snps.resize(0);
}

void SnpRepository::InitRepository(const char* file)	{
	snpCount+=CountLines(file);
	Resize(snpCount);
	snpLookup.clear();
}

//Assumes that the parser has been opened already
void SnpRepository::LoadData(GtFileParser *parser) {
	SnpAligned *snp = NULL;
	//In case this has been used once before
	parser->Reset();

	InitRepository( parser->GetIndividualCount(), parser->GetSnpCount());
	snp=parser->NextSnp();

	while (snp) {
		//cout<<"SNP acquired from parser: "<<snp->toString()<<"\n";
		if (snp->Validate()) {
			snps[curSnp++] = snp;
			snpLookup[snp->GetLabel()] = snp;

			//snps.push_back(snp);
			//curSnp++;
		}
		else	{
			pool->ReleaseSnp( snp );
			//cout<<"Not Adding "<<snp->GetID()<<" "<<snp->GetInclusionStatus()<<" \n\t\t"<<snp->toString()<<"\n";
		}
		snp=parser->NextSnp();
	}
	PostImport();
}

void SnpRepository::ParseInputFile(GtFileParser *parser, Reporting::LocusLog *locusLog) {
	SnpAligned *snp = NULL;

	//We don't know how many snps there are so far

	if (parser->Open()) {
		InitRepository( parser->GetIndividualCount(),parser->GetSnpCount());
		snp=parser->NextSnp();
		while (snp) {
			//cout<<"SNP acquired from parser: "<<snp->toString()<<"\n";
			if (snp->Validate()) {
				//cout<<"Added "<<snp->GetID()<<"\n";
				snps[curSnp++] = snp;
				snpLookup[snp->GetLabel()] = snp;

				//snps.push_back(snp);
				//curSnp++;
				if (locusLog)
					locusLog->WriteSnp(snp);
			}
			else	{
				if (locusLog)
					locusLog->WriteSnp(snp);
				pool->ReleaseSnp( snp );
				//cout<<"Not Adding "<<snp->GetID()<<" "<<snp->GetInclusionStatus()<<" \n\t\t"<<snp->toString()<<"\n";
			}
			snp=parser->NextSnp();
		}
	}
	else
		cout<<"Parser didn't open. Unable to read data\n";
	//affectedMask = *parser->GetAffectedMask();
	PostImport();
}

void SnpRepository::ParseEseBinGenofile(const char* filename) {
	SnpAligned *snp = NULL;
	BinaryGenoParser parser(filename);
	
	parser.Open();

	//We have to have someplace to put these structures. Also, we are assuming that this is being used wisely. First, initialize the 
	//Repository with ALL files then parse them
	if (snps.size() != snpCount)	{
		cout<<"* Repository not initialized. Unable to continue!\n";
		exit(1);
	}
	
	snp=parser.GetNextSnp();
	//SnpAligned *superSnp=snp;					///<Used to track down a memory problem I'm having
	while (snp) {
		snps[curSnp++] = snp;
		snpLookup[snp->GetLabel()] = snp;

		snp=parser.GetNextSnp();
	}
}

void SnpRepository::ParseBasePairTextFile(const char* filename, Reporting::LocusLog *log)	{
	//We have to have someplace to put these structures. Also, we are assuming that this is being used wisely. First, initialize the 
	//Repository with ALL files then parse them
	if (snps.size() != snpCount)	{
		cerr<<"* Repository not initialized. Unable to continue!\n";
		exit(1);
	}

	SnpAligned *snp = NULL;

	ifstream fin(filename, ios_base::in);
	char line[MAX_LINE_LENGTH];
	if (fin.eof() || fin.fail() || fin.bad() ) {
		return;
	}
	
	int individuals = 0;	
	int lineNumber=1;

	while (!fin.eof() ){
		fin.getline(line, MAX_LINE_LENGTH);
		snp=pool->GetSnp(4);
		line[strlen(line)-1]='\0';
	
		individuals=snp->ImportBpSnp(lineNumber++, curSnp, line);
		
		if (snp->Validate())	{
			log->WriteSnp(snp);
			individualCount = individuals;			//For now, let's assume this is just for genetic data. Eventually, we probably want it guideable
			snps.push_back(snp);
			curSnp++;
		}
		else {
			if (line[0] != '\0') 
				log->WriteSnp(snp);
			pool->ReleaseSnp(snp);
		}
	}
	//We don't want to misremember how far we made it
	fin.close();
}

void SnpRepository::Evaluate(SnpRecipient* repos, SnpVerificationMethod* verification) {
	SnpAligned *newSnp=NULL;
	
	for (uint curIdx=0; curIdx<snpCount; curIdx++) {
		newSnp=snps[curIdx];
		if (verification->EvaluateModel(newSnp, curIdx))
			if (repos)
				repos->Append(newSnp);
	}
}


uint SnpRepository::PerformEvaluation(SnpAligned *previousSnp, SnpRecipient *repos, SnpVerificationMethod *verification, int comboStart, int comboStop, uint spinValue) {
	uint spinner=spinValue;
	static uint totalSpins=0;
	static uint innerRuns=0;
	uint modelsSeen =0;
	SnpAligned *snp;
	
	//We have to allow for multiple searches, so allow for the static value to be reset
	if (spinValue==1) {
		totalSpins=0;
		innerRuns=0;
	}
	
	while (++spinner < snpCount) {
		snp=previousSnp->punett( snps[spinner] );
		modelsSeen++;
		//evaluate local snp for ramping...we aren't doing this yet
		if (comboStart < 1 && verification->EvaluateModel(snp, totalSpins++))
			if (repos)
				repos->Append(snp);
		//cout<<snp->GetLabel()<<" : "<<snp->GetLastMdEval()<<"\n";

		if (comboStop)
			modelsSeen+=PerformEvaluation(snp, repos, verification, comboStart - 1, comboStop-1, spinner);
		pool->ReleaseSnp(snp);
	}
	return modelsSeen;
}

uint SnpRepository::Evaluate(int comboStart, int comboEnd, SnpRecipient *repos, SnpVerificationMethod *verification) {
	SnpAligned *snp=NULL;
	uint modelsSeen = 0;
	//cout<<"Evaluation beginning\n";
	for (uint i=0; i<snpCount; i++) {
		snp=snps[i];
		modelsSeen++;
		//We can do this to always check for sub level models
		if (comboStart < 1 && verification->EvaluateModel(snp, i))
			if (repos)
				repos->Append(snp);
		//cout<<snp->GetLabel()<<" : "<<snp->GetLastMdEval()<<"\n";
		if (comboEnd>0)
			modelsSeen += PerformEvaluation(snp, repos, verification, comboStart - 1, comboEnd-1, i);
	}
	verification->PostEvaluation();
	return modelsSeen;
}

void SnpRepository::Evaluate(SnpRecipient *pass, SnpRecipient *fail, SnpVerificationMethod *verification) {
	SnpAligned *newSnp=NULL;
	
	for (uint curIdx=0; curIdx<snpCount; curIdx++) {
		newSnp=snps[curIdx];
		if (verification->EvaluateModel(newSnp, curIdx))
			pass->Append(newSnp);
		else
			fail->Append(newSnp);
	}
}

bool SnpRepository::CheckForDupes(SnpAligned *lhs, SnpAligned *rhs) {
	bool success=true;
	uint outer=lhs->GetLabelCount();
	uint inner=rhs->GetLabelCount();
	for (uint i=0; success && i<outer; i++) 
		for (uint j=0; success && j<inner; j++) 
			success=lhs->GetLabel(i) < rhs->GetLabel(j); 
	return success;
}

SnpAligned *SnpRepository::GetSingleLocusSnp(const char *snpLabel) {
	if (snpLookup.find(snpLabel) != snpLookup.end())
		return snpLookup[snpLabel];
	else {
		cout<<"Couldn't find: ."<<snpLabel<<". in the lookup ("<<snpLookup.size()<<")\n";
		map<string, SnpAligned*>::iterator itr = snpLookup.begin();
		map<string, SnpAligned*>::iterator end = snpLookup.end();

		while (itr!=end) {
			cout<<"."<<itr->first<<"."<<itr->second->GetLabel()<<"\n";

			itr++;
		}
		return NULL;
	}
}
SnpAligned *SnpRepository::GetSnp(const char *modelName) {


	char *start = (char *)modelName;					//Start of a given number in string
	char *lastChar = start + strlen(start);
	char *end = strchr(start, 'x');	//End position of the number in the string
	
	uint modelID;
	char snpID[64];
	
	SnpAligned *model = NULL;
	SnpAligned *tempModel = NULL;
	
	bool doContinue = true;
	//Build the model 1 step at a time
	while (doContinue) {
		uint len = lastChar - start;
		if (end)
			len=end-start;
		strncpy(snpID, start, len);
		snpID[len]='\0';
		modelID = atoi(snpID);
		tempModel = GetSingleLocusSnp(snpID);

		if (model) {
			tempModel = model->punett(tempModel);
			model->ReduceInstanceCount();
			model = tempModel;
		}
		else 
			model = tempModel;

		if (model)
			model->IncrementInstanceCount();
		doContinue = end != NULL;
		if (doContinue) {
			start = end + 1;
			end = strchr(start, 'x');
		}
	}

	return model;

	
}

}
