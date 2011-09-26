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
#include "sibshipfoldproduction.h"

namespace MDR  {

using namespace Genetics;
using namespace Genetics::Parser;
using namespace Genetics::Evaluation;

/**
@brief Handle the parsing of Mdr Pdt files using a buffered array to help realigning them quickly. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

class GtLineParserMdrPdt : public GtLineParser {
public:
	GtLineParserMdrPdt( int affValue, bool readLabel, FamilyRepository *famRepo, uint foldCount, int genotypeCount = -1);
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

	/**
	 * @brief Add an evaluation method to be applied to the data after loading (during post parse). 
	 * @note These will be run after load and before DSPs are constructed
	 */
	void AddPostEvaluation(FamilyRepoEvaluation *evaluator);

	uint GetIndividualCount();							///<Returns the number of individuals found
	uint GetSnpCount();									///<Returns the number of snps encountered
	void GenerateReport(ostream &os);			 		///<Generates a string report
	uint GetAffectedCount();							///<Returns the number of affected individuals
	uint GetUnaffectedCount();							///<Returns the number of unaffected individuals
	CaseControlStatus GetScrambledStatusMask();			
	void GetStatusMask(CaseControlStatus *);			///<The bit Mask produced by the affected status
	void AppendScrambledStatus(CaseControlStatus*);		///<Returns the scrambled Mask for the affected status

	bool InitData(const char *line);					///<If we decide we want to have snp labels other than line numbers

	void SetupChildren(bool doRandomize);				///<Set up any genotype data as well as build the affected status array

	CaseControlStatus *GetPgStatArray();
	uint GetPgStatArrayCount();							///<Returns the number of status objects in the parental group array
	
	void SetPedigreeStream(ostream *os);				///<Set the stream that is used for reporting on pedigree data

	void ReportOnIndividualInUse();						///<Returns the number of individuals considered during the analysis

	void SetFoldCount(uint foldCount);

	static bool WriteFoldDetails;						///<When true, details regarding the cross validation will be written to the pedigree report
	PdtFold *GetFolds(); 								///<Returns a copy of the folds (this must be deleted by the caller)

	/**
	 * @brief Marks certain pedigrees as deleted so they will be removed from the analysis
	 */
	void ExcludePedigrees(vector<string> pedigreeList);
	void ExcludeLoci(vector<string> lociList);

	/**
	 * @brief Write the repository to file
	 */
	void WritePedFile(const char *filename);
	
	/**
	 * @brief Writes a copy of the dataset but only the the loci specified by modelLoci (format 5x10)
	 */
	void WriteModelLociOnly(const char *modelLoci, const char *filename);
	void WriteSaScript(const char *filename, vector<SnpAligned *>, const char *targetPed, bool withInteraction);
protected:
	vector<FamilyRepoEvaluation *> postEvalMethods;		///<These are used to process the family data prior to setting up DSPs

	FamilyRepository *famRepo;							///<The source for the genetic data (and family information)
	vector<GenotypeData *> includedChildren;			///<The genotype data to be used for the snps
	SnpPool *pool;										///<Snp pool used to acquire and release items
	vector<SnpAligned::Label> labels;					///<The array of labels
	int affectedValue;									///<If our data has unknown status, this has to be changed - This seems to be acceptable, but users should be aware to clean their data before hand
	bool readLabel;										///<Indicates if there is expected to be a header
	int unknownGenotype;								///<This is the value used to assign
	uint totalIndividualsSeen;							///<The number of individuals that have been passed to the parser from the file
	CaseControlStatus individualsInUse;					///<Individuals in use (not duplicated)
	uint pgStatArrayCount;

 	CaseControlStatus overAllStatus;					///<This is the status for the individuals regardless of the family structure
	ostream *pedStream;									///<Used to write to the pedigree log

	CaseControlStatus *pgStatArray;
	uint foldCount;										///<The number of cross validation folds
	PdtFold *folds;										///<Array of folds containing the necessary pieces for evaluation
	vector<string>	pedExclusionList;					///<List of pedigrees that are to be ignored
	vector<string> lociExclusionList;
	BitSetType snpExclusions;							///<Mask indicating which snps should be ignored
	uint lociObserved;
};


inline
void GtLineParserMdrPdt::SetFoldCount(uint folds) {
	foldCount = folds;
}

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
uint GtLineParserMdrPdt::GetPgStatArrayCount() {
	return pgStatArrayCount;
}


inline
CaseControlStatus *GtLineParserMdrPdt::GetPgStatArray() {
	return pgStatArray;
}


inline
PdtFold *GtLineParserMdrPdt::GetFolds() {
	PdtFold *newFolds = new PdtFold[foldCount];
	for (uint i=0;i<foldCount; i++)
		newFolds[i]=folds[i];
	return newFolds;
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
GtLineParserMdrPdt::GtLineParserMdrPdt( int affValue, bool readLabel, FamilyRepository *famRepo, uint foldCount, int genotypeCount/*=-1*/) :
			GtLineParser(genotypeCount), famRepo(famRepo), affectedValue(affValue), 
			readLabel(readLabel), totalIndividualsSeen(0), pgStatArrayCount(0), 
			pedStream(NULL), pgStatArray(NULL), foldCount(foldCount), folds(NULL), lociObserved(0) {
	pool=SnpPool::Instance();
}

inline
GtLineParserMdrPdt::~GtLineParserMdrPdt() {
	if (pgStatArray)
		delete[] pgStatArray;

	if (folds)
		delete[] folds;
}


}


#endif
