//
// C++ Interface: snprepostxtsorted
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESESNPREPOSTXTSORTED_H
#define ESESNPREPOSTXTSORTED_H
#include "snprecipient.h"

#include "snpcontainer.h"
#include "permutationtest.h"

namespace Genetics {


/**
@brief Class to use for building basic text reports

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpReposTxtSorted : public SnpRecipient, public SnpReporter
{
public:
struct ReportModel {
	ReportModel() : labelCount(0), score(0) {}

	ReportModel(SnpAligned *snp) {
		score=snp->GetLastMdEval();
		labelCount=snp->GetLabelCount();
		
		for (uint i=0; i<labelCount; i++) 
			labelArray.push_back(SnpAligned::Label(snp->GetLineNumber(i), snp->GetLabel(i)));
	}

	ReportModel(const ReportModel& other) {
		labelCount=other.labelCount;
		score=other.score;
		
		for (uint i=0; i<labelCount; i++) 
			labelArray.push_back(other.labelArray[i]);
	}

	string GetLabel() {
		char lineNumber[2048];
		sprintf(lineNumber, "%d", labelArray[0].label);
	
		for (uint i=1; i<labelCount; i++)
			sprintf(lineNumber, "%sx%d", lineNumber, labelArray[i].label);

		return lineNumber;
	}
	

	string GetLineNumber() {
		char lineNumber[2048];
		sprintf(lineNumber, "%d", labelArray[0].lineNumber);
	
		for (uint i=1; i<labelCount; i++)
			sprintf(lineNumber, "%sx%d", lineNumber, labelArray[i].lineNumber);

		return lineNumber;
	}

	
	vector<SnpAligned::Label> labelArray;
	uint labelCount;
	float score;

	bool operator<(const ReportModel& other) {
		return score<other.score;
	}

	~ReportModel() {
	}
};
struct ReportModelLT {
	bool operator()(const ReportModel& s1, const ReportModel& s2) {
		return s2.score < s1.score;
	}
};

struct ReportEntry {
	string shortEntry;
	string longEntry;
	string id;
	uint val;
	ReportEntry(string id, uint val, string s, string l) :  shortEntry(s), longEntry(l), id(id), val(val) {}
	bool operator<(const ReportEntry& other) { 
		return val<other.val;
	}
};
struct ReportEntryLT {
	bool operator()(const ReportEntry& s1, const ReportEntry& s2) {
		return s2.val < s1.val ;
	}	
	bool operator()(const ReportModel& s1, const ReportModel& s2) {
		return s2.score < s1.score;
	}

};

typedef vector <ReportModel> ReportArray;



 	/**
	 * @brief Construct the container
	 * @param snpCount indicates how many snps are to be able to be stored
	 */
    SnpReposTxtSorted(uint snpCount);

    ~SnpReposTxtSorted();
	/**
	 * @brief Sticks a new model into the array. 
	 * The array will grow by a predefined amount in the event of adding more snps than had been initially
	 * allowed for
	 */
	void Append(SnpAligned *newSnp);

	/**
	 * @brief Get the number of snps inside the repository
	 */
	uint GetEntryCount(uint loci);

	uint GetTotalEntryCount(uint loci);

	/**
	 * @brief Return the snp at position n
	 */
	string GetReportEntry(uint loci, uint n);
	ReportModel *GetReportModel(uint loci, uint n){
		if (n<contents[loci]->size() == 0)  
			return NULL;
	
		ReportModel &model = contents[loci]->at(n);
		return &model;
	}

	/**
	 * @brief Resize the local content array.
	 */
	void ResizeArray(uint loci);

	/**
	 * @brief perform a sort on each of the sub report arrays
	 */
	void Sort();

	Reporting::PermutationTestDist *distribution;		///<The distribution of permutation tests

protected:
	vector< ReportArray *> contents;
	uint snpCount;										///<The number of snps to be stored
};

inline
SnpReposTxtSorted::SnpReposTxtSorted(uint snpCount) : distribution(NULL), snpCount(snpCount) {
}

inline
SnpReposTxtSorted::~SnpReposTxtSorted() {
	uint count=contents.size();
	ReportArray *reports;
	for (uint i=0; i<count; i++)	{
		reports=contents.back();
		contents.pop_back();
		delete reports;
	}
	
}

inline
uint SnpReposTxtSorted::GetTotalEntryCount(uint loci) {
	if (contents.size() > loci)
		return contents[loci]->size();
	else
		return 0;
}

inline
uint SnpReposTxtSorted::GetEntryCount(uint loci) {
	if (contents.size() < 1)
		return 0;

	uint count=contents[loci]->size();

	if (snpCount == 0 || count < snpCount)
		return count;
	else
		return snpCount;
}

/**
 * @todo SnpReposTxtSorted::GetReportEntry() should return an entry and the report should contruct entries...
 */
inline
string SnpReposTxtSorted::GetReportEntry(uint loci, uint snpIdx) {
	if (snpIdx<contents[loci]->size() == 0)  
		return "";

	char basicEntry[128];
	ReportModel &model = contents[loci]->at(snpIdx);
	double pVal=0.0;
	if (distribution)
		pVal=distribution->GetPValue(model.labelCount, model.score);
	sprintf(basicEntry, "%s\t%.4f %.4f", model.GetLabel().c_str(), model.score, pVal);
	return basicEntry;
}

inline
void SnpReposTxtSorted::Sort() {
	uint count=contents.size();
	for (uint i=0; i<count; i++)
		sort(contents[i]->begin(), contents[i]->end(),ReportEntryLT());
}

inline
void SnpReposTxtSorted::ResizeArray(uint loci) {
	//ReportArray *reports;
	for (uint i=contents.size(); i< loci; i++) 
		contents.push_back( new ReportArray() );
}

}

#endif
