//
// C++ Interface: snprepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef E_SESNPREPOSITORY_H
#define E_SESNPREPOSITORY_H

#include <vector>
#include "snppool.h"
#include "snprecipient.h"
#include "snpverificationmethod.h"
#include "gtfileparser.h"
#include "locuslog.h"

namespace Genetics {

using namespace Genetics;
using namespace Parser;

using namespace ValidationFunctors;

typedef vector<SnpAligned *> SnpArray;
#define DEFAULT_GROWTH_SIZE 150
/**
This contains and provides the storage for the various snps

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpRepository : public SnpRecipient
{
public:
    SnpRepository();

    ~SnpRepository();

	/**
	 * @brief Just runs through the file and counts the number of lines
	 */
	int CountLines(const char *file);

	/**
	 * @brief Return the number of individuals have been encountered. 
	 * This will either reflect the last number of individuals encountered or the uniform individual count
	 * after a validation has been performed. 
	 */
	uint GetIndividualCount();

	 /**
	  * Returns the number of snps associated the repository
	  */
	uint GetSnpCount();

	/**
	 * @brief Clears the idividual/snp counts;
	 */
	void InitRepository(uint indCount, uint snCount);

	/**
	 * @brief Parses the file to estimate the snp counts found inside the file
	 * @note At this time, the indivuals aren't counted. It will be a validation step that is used to 
	 * recognize that some snps aren't all sized the same
	 */
	void InitRepository(const char* file);

	/**
	 * @brief Loads in the genotypes from a base pair encoded text file and stores them inside the SnpArray
	 * @param file File that is to be loaded
	 */
	void ParseBasePairTextFile(const char* file, Reporting::LocusLog *locusLog);

	/**
	 * @brief Parses an ese binary file
	 */
	void ParseEseBinGenofile(const char* filename);

	/**
	 * @brief Uses parser to populate the snp data from the file
	 * @param parser This is the helper object that knows how to parse the file and populate the repository with snps
	 */
	void ParseInputFile(GtFileParser *parser, Reporting::LocusLog *locusLog);
	void LoadData(GtFileParser *parser);

	/**
	 * @brief Called when all in put files have been read, this sets the snpCount appropriately
	 */
	void PostImport();

	
	/**
	 * @brief Evaluates each element in the current repository and adds a copy to the new one if it passes validation
	 * @param repos New repository where the succesful candidates go
	 * @param verification Functor used to evaluate the local snps
	 */
	void Evaluate(SnpRecipient* repos, SnpVerificationMethod* verification);

	/**
	 * @brief Splits the repository into two different repositories based on the outcome of verification
	 * @param pass Recieves each of the snps that passes verification
	 * @param fail Recieves each of the snps that fails the verification
	 * @param verification Functor used to evaluate local snps
	 */
	void Evaluate(SnpRecipient *pass, SnpRecipient *fail, SnpVerificationMethod *verification);
	uint PerformEvaluation(SnpAligned *previousSnp, SnpRecipient *repos, SnpVerificationMethod *verification, int comboStart, int comboEnd, uint spinValue);
	uint Evaluate(int combStart, int comboEnd, SnpRecipient *repos, SnpVerificationMethod *verification);
	/**
	 * @brief Perform a search for the perfect model in this repository
	 */
	//void SearchForPerfectModels(SnpRecipient* newRepo);	 
	/**
	 * Returns the SNP associated with the <I>idx</I>th position in the repository
	 */
	SnpAligned *GetSnp(uint idx);

	/**
	 * @brief Sticks a new model into the array. 
	 * The array will grow by a predefined amount in the event of adding more snps than had been initially
	 * allowed for
	 */
	void Append(SnpAligned *newSnp);

	/**
	 * @brief Sets the amount of snps that added to the array each time we append 1 more than we have space for
	 */
	void SetGrowby(uint amount);

	//BitSetType *GetAffMask();

	SnpAligned *GetSnp(const char *modelID);

	SnpAligned *BuildRandomModel(uint modelSize);
	SnpAligned *GetSingleLocusSnp(const char *modelID);
	void PurgeSnps();
	//void Cross(int startIndex, SnpAligned *snp, SnpRecipient* repos, SnpVerificationMethod& verification);
protected:

	void Grow(int count);
	
	/**
	 * @brief change the size of the structure putting NULLs into the new positions and releasing those that are in excess of the size. 
	 * (The ones released are those at the end of the array)
	 * @param count The new size of the array
	 */
	void Resize(int count);

	bool CheckForDupes(SnpAligned *lhs, SnpAligned *rhs);

	uint curSnp;						///<The current snp index. Used during loading
	uint individualCount;				///<The number of individuals the snps are thought to contain. This is only reliable after validation has been performed
	uint snpCount;						///<The number of snps in the repository (or the number that is expected to be found)
	SnpArray snps;						///<The actual snps
	SnpPool *pool;						///<Used to maintain the pool of available snps
	uint growby;						///<Used to increase the size of the array during appending
	map<string, SnpAligned*> snpLookup;	///<This is used by GetSnp(char *) to return the correct snp
	//BitSetType affectedMask; 			///<This is used to represent the cases (1s)
	//BitSetType ignoreMask;				///<This is used to store which individuals are not to be used in analysis
};

inline
void SnpRepository::PostImport() {
	snpCount=curSnp;
}

inline 
void SnpRepository::SetGrowby(uint amount) {
	growby=amount;
}

inline
void SnpRepository::Append(SnpAligned *newSnp) {
	if (snpCount == 0 || snpCount > snps.size() - 1)
		Resize(snpCount+growby);
	newSnp->IncrementInstanceCount();
	//cout<<"Appending Snp: "<<newSnp->GetLabel()<<"\n";
	//cout<<"Appending snp("<<snpCount<<" <- "<<newSnp->GetID()<<"\n";
	snps[snpCount++] = newSnp;
	snpLookup[newSnp->GetLabel()] = newSnp;
	curSnp=snpCount;
}

inline
void SnpRepository::InitRepository(uint indCount, uint snCount)	{
	curSnp=0;
	PurgeSnps();

	individualCount=indCount;
	snpCount=snCount;
	Resize(snpCount);

}

inline
void SnpRepository::Grow(int count) {
	int curSize=snps.size();
	if (curSize < count -1)
		Resize(curSize+growby);
}
/*
inline
BitSetType *SnpRepository::GetAffMask() {
	return &affectedMask;
}
*/

inline
SnpAligned *SnpRepository::GetSnp(uint idx) {
	if (idx<snpCount)
		return snps[idx];
	else 	
		return NULL;
}

inline
SnpAligned *SnpRepository::BuildRandomModel(uint modelSize) {
	uint lociCount = GetSnpCount();
	BitSetType rndLoci(lociCount);
	stringstream modelID;

	uint n = Utility::Random::globalGenerator((int)lociCount);
	modelID<<n + 1;
	for (uint m = 1; m<modelSize; m++) {
		rndLoci[n] = true;
		while (rndLoci[n]) 
			n = Utility::Random::globalGenerator((int)lociCount);
		modelID<<"x"<<n + 1;
	}
	return GetSnp(modelID.str().c_str());
}

inline
SnpRepository::SnpRepository()	: curSnp(0), individualCount(0), snpCount(0), growby(DEFAULT_GROWTH_SIZE)
{
	pool=SnpPool::Instance();
}

inline
uint SnpRepository::GetIndividualCount()	{
	return individualCount;
}

inline
uint SnpRepository::GetSnpCount()	{
	return snpCount;
}


}

#endif
