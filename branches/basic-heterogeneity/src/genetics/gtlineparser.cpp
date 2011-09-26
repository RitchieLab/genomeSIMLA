//
// C++ Implementation: gtlineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gtlineparser.h"
#include "utility/utility.h"


namespace Genetics {

namespace Parser {


void GtLineParserBP::GenerateReport(ostream &os) {
	os<<"\tInput File Type: Base Pair Format\n";
	os<<"\tIndividuals (aff: "<<statusMask.affected.count()<<", unaff: "<<statusMask.unaffected.count()<<")\n";
	os<<"\tAffected Mask: "<<statusMask.affected<<"\tUnaffected Mask: "<<statusMask.unaffected<<"\n";
}

/**
 * Since MDR uses 0, 1, 2 for genotypes, we need to adjust our values up 1. 
 */
int GtLineParserMDR::ConvertData(const char *value) {
	int v = atoi(value);
	if (v==unknownGenotype)
		return 0;

	return v + 1;
}

void GtLineParserBP::AppendScrambledStatus(CaseControlStatus *status) {
	*status=statusMask.MakeRandomCopy();
}

void GtLineParserBP::GetStatusMask(CaseControlStatus *status) {
	*status=statusMask;
}

bool GtLineParserMDR::ParseLine( const char *line, uint val) {

	if (val==0) {
		if (InitData(line))
			return false;
	}

	if (strlen(line) < 1) { return false; }

	vector<string> tokens;
	Utility::TokenizeString( line, tokens, " \t,");
	
	vector<string>::iterator i=tokens.begin();
	vector<string>::iterator end = tokens.end();

	int snCount =0;
	string newValue;
	vector<int> indData;
	if (i != end) {
		//First we need to grab the affected status
		indData.push_back(atoi((i++)->c_str()));
		//Iterate over the remaining snps bits
		for (; i != end; ++i) {	
			snCount++;
			newValue=*i;
			//indData.push_back(1);
			indData.push_back(ConvertData(newValue.c_str()));
		}
		
	}
	if (indData.size() > 0)	{
		data.push_back(indData);
		return true;
	}
	else
		return false;
}
/*
bool GtLineParserMDR::ParseLine( const char *line, uint val) {
	//If this is the first column, we need to populate our members
	if (val==0)		{
		if (InitData(line)) {
			return false;
		}
	}
	
	if (strlen(line) < 1) {
		return false;
	}
	
	string convertedLine=line;
	boost::char_separator<char> sep(" \t,", "", boost::drop_empty_tokens);
	strtokenizer tok(convertedLine, sep);						
	strtokenizer::iterator i = tok.begin();
	strtokenizer::iterator end = tok.end();

	int snCount =0;
	string newValue;
	vector<int> indData;
	if (i != end) {
		//First we need to grab the affected status
		indData.push_back(atoi((i++)->c_str()));
		//Iterate over the remaining snps bits
		for (; i != end; ++i) {	
			snCount++;
			newValue=*i;
			//indData.push_back(1);
			indData.push_back(ConvertData(newValue.c_str()));
		}
		
	}
	if (indData.size() > 0)	{
		data.push_back(indData);
		return true;
	}
	else
		return false;
}
*/
void GtLineParserMDR::GenerateReport(ostream &os) {
	os<<"\tInput File Type: MDR Format\n";
	//CaseControlStatus st=statusMask->CombinedStatus();
	
	os<<"\tIndividuals (aff: "<<stMask.affected.count()<<", unaff: "<<stMask.unaffected.count()<<")\n";
	os<<"\tAffected Mask: "<<stMask.affected<<"\tUnaffected Mask: "<<stMask.unaffected<<"\n";
}




bool GtLineParserMDR::InitData(const char *line) {
	vector<string> tokens;
	Utility::TokenizeString( line, tokens, " \t,");
	
	vector<string>::iterator i=tokens.begin();
	vector<string>::iterator end = tokens.end();
	
	int snpCount=-1;
	
	
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

void GtLineParserMDR::GetStatusMask(CaseControlStatus *status) {
	uint indCount=data.size();
	BitSetType aff(indCount);
	BitSetType unaff(indCount);
	
	//Is there a way for an MDR formatted file to have unknown status?
 	//People are saying that there isn't. But, there is the possibility for us 
	//to chop this status up by slicing for cross validation
	for (uint indID=0; indID<indCount; indID++)
		aff[indID]=data[indID][0]==affectedValue;
	unaff = ~aff;
	
	stMask = CaseControlStatus(aff, unaff);
	*status=stMask;
}

void GtLineParserMDR::AppendScrambledStatus(CaseControlStatus *status) {
	//MDR is a very simple method- just acquire a single randomized copy of the mask
	*status=stMask.MakeRandomCopy();
}




void GtLineParserMDR::PrintData() {
	for (uint x=0; x<data.size(); x++) {
		for (uint y=0; y<data[x].size(); y++)
			cout<<data[x][y]<<".";
		cout<<"\n";
	}
	cout<<"\n";
}

SnpAligned *GtLineParserMDR::GetSnp(uint id) {
	uint snpID=id+1;
	SnpAligned *snp=NULL;
	assert(data[0].size() == data[1].size()); 
	if (snpID < data[0].size()) {
		uint indCount=data.size();
		snp=pool->GetSnp(genoCount, indCount);
		snp->Reset();
		//snp->SetIndividualCount(indCount);
		for (uint indID=0; indID<indCount; indID++)	{
			assert(data[indID].size() > snpID);
			snp->SetValue(indID, data[indID][snpID]);
		}
		snp->SetLabel(labels[snpID]);
		BitSetType &missing = snp->GetGenotype( 0 );
		missing=~missing; 
	}
	return snp;			
}


}

}
