//
// C++ Interface: individual
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTINDIVIDUAL_H
#define MDRPDTINDIVIDUAL_H
#include <string>
#include <vector>

namespace MdrPDT {

/**
	Information regarding a family member in a pedigree dataset. This is used to 
	cache the genotye values, since we will need to rewrite them for each PTest 
	(or to revert back to the real evaluation)

	Construction/Deletion		Individuals are created/destroyed solely by the 
	pedigree, since the pedigree must be responsible for mantaining their association. 

	Affected status is used only for the real status. During permutation, it is expected
	that the object writing the genotypes to the repository are handling the status 
	properly
	
	Usage:
		(From Pedigree)		Instantiate, Set Pedigree details, Initialize virtual siblings
		(everywhere else)	Get genotypes based on index and genotype counts. -1 is missing as 
							well as 

	@author 
*/
class Individual{
friend class Pedigree;
public:

	/**
	 * @brief Returns t/f regarding affected status
	 */
	bool IsAffected();

	/**
	 * @brief sets the affected status for the user
	 */
	void IsAffected(bool isAffected);

	/**
	 * @brief Adds the genotype to the member's list
	 */
	void AddGenotype(char val);

	/**
	 * @brief Gets the genotype at a given index
	 */
	char GetGenotype(size_t idx);

	void SetGenotype(size_t idx, char gt);

	/**
	 * @brief returns the genotype for the non-transmitted sibling at index, idx
	 */
	char GetVirtualGenotype(size_t idx);

	/**
	 * @brief Returns the IDs for the parents
	 */
	void GetParentIDs(std::string& pat, std::string& mat);

	Individual *GetMother() { return mother; }
	Individual *GetFather() { return father; }

	void SetStatus(char status);
	
	void SetMother(Individual*m) {mother=m;}
	void SetFather(Individual*f) {father=f;}

	/**
	 * @brief We don't care about this value, but we'll hold onto it
	 */
	void SetGender(char gender);
	/**
	 * @brief sets the IDs for the parents
	 */
	void SetParentIDs(const std::string& pat, const std::string& mat);	

	/**
	 * @brief Returns the ID for the member
	 */
	std::string GetID();
	
	/**
	 * @brief Sets the ID for the member
	 */
	void SetID(const char *id);

	/**
	 * @brief returns the number of genotypes for this individual
	 */
	size_t CountGenotypes();

	/**
	 * @brief Sets up the nontransmitted allele
	 */
	bool SetupVirtuals(std::ostream& os);

	/**
	 * @brief Indicate if there are virtual genotypes or not
	 */
	bool HasVirtualGentoypes();

	/**
	 * @brief Returns the integer pedigree ID. 
	 * This is guaranteed to be an integer and is no bigger than the number of 
	 * individual pedigrees observed in the dataset. This value is probably
	 * different from the one listed in the dataset itself and should not 
	 * be reported to the user
	 */
	int GetPedigreeID();

	void GenerateReport(std::ostream& os, bool reportVirtuals);

	void ReportGenotypes(std::ostream& os, bool reportVirtual, std::string sep);
	/**
	 * @brief write the genotypes in use to stream. 
	 * @param os The stream to be written to
	 * @param writeVirtuals If true, the virtual sibling will be written instead
	 */
	void Write(std::ostream& os, bool writeVirtuals=false);
	/**
	 * @brief Indicate that an individual should be ignored during analysis
	 */
	void DropFromAnalysis(bool doDrop) { dropFromAnalysis=doDrop; }
	/**
	 * @brief Returns true if the individual is not part of anlaysis
	 */
	bool DropFromAnalysis() { return dropFromAnalysis; }

	static char AffectedValue;
	static char UnaffectedValue;
	static bool UseMerlin;
	int GetValidGenotypeCount() { return validGenotypes; }
protected:
	/**
	 * @brief All initialization is done by the pedigree
	 */
    Individual(const char *id, const char *pedID, int pedigreeIdx);

    ~Individual();

	bool dropFromAnalysis;



	std::string id;						///<Individual ID
	std::string pedID;					///<Pedigree ID
	std::string pat;					///<Paternal ID
	std::string mat;					///<Maternal ID

	Individual *father;					///<Pointer to the father object (if exsts)
	Individual *mother;					///<Pointer to the mother
	int pedigreeID;						///<Index of the pedigree within all others in the dataset

	int status;							///<Affected Status (0 - unaffected, 1 - affected)
	char gender;						///<Gender, in case we need to write it back to file
	std::vector<char > genotypes;		///<Genotype information in (0, 1, 2)
	std::vector<char > virtualGenotypes;	///<Non-Transmitted alleles, based on parents. This is only populated if the parents exist in the dataset

	int validGenotypes;					///<Number of present genotype values
};

}

#endif
