//
// C++ Interface: gtlineparsermdrpdt
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_PARSERGTLINEPARSERMDRPDT_H
#define GENETICS_PARSERGTLINEPARSERMDRPDT_H

#include "genetics/gtlineparser.h"
#include "genetics/familyrepository.h"
#include "genetics/familymember.h"
#include "genetics/familyevaluation.h"

namespace ESE  {

using namespace Genetics;
using namespace Genetics::Parser;
using namespace Genetics::Evaluation;

/**
@brief Handle the parsing of Mdr Pdt files using a buffered array to help realigning them quickly. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

class GtLineParserMdrPdt : public GtLineParser {
public:
	GtLineParserMdrPdt( int affValue, bool readLabel, FamilyRepository *famRepo, int genotypeCount = -1);
	~GtLineParserMdrPdt();

	/**
	 * @brief called by the file parser object with a single line and the line number
	 */
	bool ParseLine(const char *line, uint val);

	/**
	 * @brief Returns the next snp
	 */
	SnpAligned *GetSnp(uint snpID);
	
	/**
	 * @brief this is used to prep the incudedChildren vector
	 * This is called after the parser has completed all of the input files
	 */
	void PostParse();

	
	uint GetIndividualCount();							///<Returns the number of individuals found
	uint GetSnpCount();									///<Returns the number of snps encountered
	void GenerateReport(ostream &os);			 		///<Generates a string report
	uint GetAffectedCount();							///<Returns the number of affected individuals
	uint GetUnaffectedCount();							///<Returns the number of unaffected individuals
	//CaseControlStatus GetStatusMask();					///<The bit Mask produced by the affected status
	CaseControlStatus GetScrambledStatusMask();			
	void GetStatusMask(CaseControlStatus *);			///<The bit Mask produced by the affected status
	void AppendScrambledStatus(CaseControlStatus*);		///<Returns the scrambled Mask for the affected status

	bool InitData(const char *line);					///<If we decide we want to have snp labels other than line numbers

	void SetupChildren(bool doRandomize);				///<Set up any genotype data as well as build the affected status array

	CaseControlStatus *GetStatArray();					///<Return the status array
	CaseControlStatus *GetPgStatArray();
	uint GetStatArrayCount();							///<Return the number of status objects in the array
	uint GetPgStatArrayCount();							///<Returns the number of status objects in the parental group array
	
	vector<CaseControlStatus> *GetStatusArray();		///<Returns the set of masks associated with the various families. Valid only after "GetSnp(idx) has been called
	void SetPedigreeStream(ostream *os);				///<Set the stream that is used for reporting on pedigree data

	void ReportOnIndividualInUse();						///<Returns the number of individuals considered during the analysis
	void DoCreateVirtualSibs(bool doCreateSibs);		///<Indicates that the user does want to create virtual siblings
protected:
	FamilyRepository *famRepo;							///<The source for the genetic data (and family information)
	vector<GenotypeData *> includedChildren;			///<The genotype data to be used for the snps
	SnpPool *pool;										///<Snp pool used to acquire and release items
	vector<SnpAligned::Label> labels;					///<The array of labels
	int affectedValue;									///<If our data has unknown status, this has to be changed - This seems to be acceptable, but users should be aware to clean their data before hand
	bool readLabel;										///<Indicates if there is expected to be a header
	int unknownGenotype;								///<This is the value used to assign
	uint totalIndividualsSeen;							///<The number of individuals that have been passed to the parser from the file
	CaseControlStatus individualsInUse;					///<Individuals in use (not duplicated)
	uint statArrayCount;								///<This is the size of the status array
	uint pgStatArrayCount;
	vector<CaseControlStatus> statusArray;				///<Used to keep up with the various case control sets overall
 	CaseControlStatus overAllStatus;					///<This is the status for the individuals regardless of the family structure
	ostream *pedStream;									///<Used to write to the pedigree log
	CaseControlStatus *statArray;						///<Array of statuses (to replace the vector)
	CaseControlStatus *pgStatArray;
	bool doCreateVirtualSibs;							///<Indicates that we do or do not expand triads
};


inline
uint GtLineParserMdrPdt::GetIndividualCount() {
	return includedChildren.size();
}

inline
uint GtLineParserMdrPdt::GetSnpCount() {
	if (includedChildren.size() > 0)
		return includedChildren[0]->GetGenotypeNumber(0);
	else
		return 0;
}

inline
void GtLineParserMdrPdt::SetPedigreeStream(ostream *os) {
	pedStream=os;
}

inline
uint GtLineParserMdrPdt::GetStatArrayCount() {
	return statArrayCount;
}

inline
uint GtLineParserMdrPdt::GetPgStatArrayCount() {
	return pgStatArrayCount;
}


inline
CaseControlStatus *GtLineParserMdrPdt::GetStatArray() {
	return statArray;
}

inline
CaseControlStatus *GtLineParserMdrPdt::GetPgStatArray() {
	return pgStatArray;
}


inline
vector<CaseControlStatus> *GtLineParserMdrPdt::GetStatusArray() {
	return &statusArray;
}


inline
uint GtLineParserMdrPdt::GetAffectedCount() {
	CaseControlStatus status;
	GetStatusMask(&status);
	return status.affected.count();
}

inline
uint GtLineParserMdrPdt::GetUnaffectedCount() {
	CaseControlStatus status;
	GetStatusMask(&status);
	return status.unaffected.count();
}

inline
void GtLineParserMdrPdt::DoCreateVirtualSibs(bool doCreate) {
	doCreateVirtualSibs=doCreate;
}


inline
void GtLineParserMdrPdt::GetStatusMask(CaseControlStatus *status) {
	*status = overAllStatus;
}



inline
void GtLineParserMdrPdt::AppendScrambledStatus(CaseControlStatus *status) {
	//Iterate over the families
 	cout<<"Scrambled status doesn't work for MdrPDT parser. Use SetUpChildren(true)!\n"; 
	assert(0);

}



inline
void GtLineParserMdrPdt::PostParse() {
	SetupChildren(false);
}






inline
GtLineParserMdrPdt::GtLineParserMdrPdt( int affValue, bool readLabel, FamilyRepository *famRepo, int genotypeCount/*=-1*/) :
			GtLineParser(genotypeCount), famRepo(famRepo), affectedValue(affValue), 
			readLabel(readLabel), totalIndividualsSeen(0), statArrayCount(0), 
			pgStatArrayCount(0), pedStream(NULL), statArray(NULL), pgStatArray(NULL), doCreateVirtualSibs(false) {
	pool=SnpPool::Instance();
}

inline
GtLineParserMdrPdt::~GtLineParserMdrPdt() {
	delete[] statArray;
	if (pgStatArray)
		delete[] pgStatArray;
}


}


#endif
