//
// C++ Interface: genotyperepository
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTGENOTYPEREPOSITORY_H
#define MDRPDTGENOTYPEREPOSITORY_H
#include <iostream>
#include <vector>
#include "utility/types.h"

namespace MdrPDT {

/**
 * @brief This is just a place to hold the genotype information and help obtain the corect pointers
 * It is important that all users parse the row properly. The repository 
	just releases a pointer and lets the client handle it however it wishes...which
	means the client can really get messed up, if it is not synched up with 
	the layout of the data
 */
class GenotypeRepository {
public:
	GenotypeRepository();	
	GenotypeRepository(const GenotypeRepository& other);
	~GenotypeRepository();

	/**
	 * @brief Initialize the repository with the correct number of SNPs and Individuals
	 */
	char *Initialize(int snpCount, int indCount);
	/**
 	 * Helper function which will return the right spot within the genotypes 
	 * array for a given SNP (assuming that all snps in the index are present)
     */
	char *GetSNP(int idx);
	char *GetSNPs() { return genotypes; }
	/**
	 * @brief Returns the number of snps
	 */
	int GetSnpCount() { return snpCount; }

	/**
	 * @brief Returns the number of individuals
	 */
	int GetIndividualCount() { return indCount; }

	void Dump(std::ostream& os);
	char *GetFolds();
	int *GetPedigreeIDs();

	void InitExclusionList(std::vector<std::string>& list);

	std::string GetGenotypeEncoding(int idx, int genotype);

	void DoAnalyzeSNP(int idx, bool doAnalyze);
	bool DoAnalyzeSNP(int idx);
	std::vector<std::string>* GetSnpLabels() { return snpLabels;}
	void SetSnpLabels(std::vector<std::string>* labels) { snpLabels=labels; }
	std::vector<int>* GetDSPCounts(){return &dspCount; }
	void SetDSPCount(int fold, int count);
	void SetFoldCount(int count);
	int GetFoldCount();

	/**
	 * @brief Counts missing for each fold and reports them, marking loci which exceed threadhold to be dropped from the analysis.
	 */
	void EvaluateMissingData(float maxMissingThreshold, std::ostream& report);
protected:
	char *genotypes;						///<Real Genotypes 
	int *pedIDs;							///<Extra stuff, like pedigree IDs
	char *folds;							///<Folds

	/**
	 * This is the width of a row. This might be somewhat larger than the actual 
	 * number of individuals, due to expanding trios and whatnot. The calling routine
	 * should take this into account (also, we will add extra values for headers, if
	 * we need them)
 	 */ 
	int rowWidth;
	int indCount;							///<Number of individuals in the dataset
	int snpCount;							///<Number of SNPs

	/**
	 * @brief T/F indicating that a snp is NOT on the ecluded list
	 */
	Utility::BitSetType doAnalyzeSNP;
	std::vector<int> dspCount;				///<Number of DSPs for each fold
	std::vector<std::string>* snpLabels;	///<Labels associated with each of the loci
};

inline
bool GenotypeRepository::DoAnalyzeSNP(int idx) {
	if (idx>0 && idx <=(int)doAnalyzeSNP.size())
		return doAnalyzeSNP[idx-1];
	return false;
}

inline
void GenotypeRepository::DoAnalyzeSNP(int idx, bool doAnalyze) {
	doAnalyzeSNP[idx-1]=doAnalyze;
}
}

#endif
