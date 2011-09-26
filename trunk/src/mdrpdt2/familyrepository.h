//
// C++ Interface: familyrepository
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTFAMILYREPOSITORY_H
#define MDRPDTFAMILYREPOSITORY_H
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include "utility/lineparser.h"
#include "pedigree.h"
#include "utility/random.h"
#include "genotyperepository.h"
#include <iostream>
namespace MdrPDT {


/**
@brief Simple conversion method

	@author Eric Torstenson
*/

//Later, we can expand this using something more interesting by using inheritance
class GenotypeConversion {
public:
	int operator()(const char *al1, const char *al2, int idx) {
		int genotype = (atoi(al1) - 1) + (atoi(al2) - 1) + 1;
		if (atoi(al1) < -1 or atoi(al1) > 2 || atoi(al2) < -1 or atoi(al2) > 2) {
			std::cerr<<"Invalid genotypes have been found: "<<al1<<" "<<al2<<"\n";
			std::cerr<<"\tGenotypes should be encoded as 1 1, 1 2, or 2 2 (missing as 0 0, or -1 -1)\n";
			assert(atoi(al1) > -2 and atoi(al1) < 3 && atoi(al2) > -2 or atoi(al2) < 3);
		}
		return genotype;
	}
};

/**
 * @brief Help to keep up with how many each pedigree will contribute
 */
struct FamilyRecord {
	int offset;				///<This is the index of the first member in the raw data
	int dspCount;			///<The number of pairs to be written to the array
	std::string pedigreeID;		///<Useful for tracking down who is in which fold
	
	FamilyRecord(int offset, int dspCount, const char *pedID) : offset(offset), dspCount(dspCount), pedigreeID(pedID) { }
	
	bool operator<(const FamilyRecord& other) const {
		return dspCount > other.dspCount;
	}
};
/**
	@author Eric Torstenson
	@brief Manages the pedigree information, as well as structures the genotype data 
	Usage:
		For each dataset to be analyzed:
			Load file
			Initialize Cross Validation (defaults to 1 fold if this isn't called)
			InitializeData (returns the raw data to be analyzed
	 		For each PTest, InitPTest(seed) (returns the raw data to be analyzed)
*/

class PedigreeRepository : public Utility::AsciiParser{
public:
	PedigreeRepository();
    ~PedigreeRepository();


	/**
	 * @brief This allows basic iteration over the contents of the repository. 
	 * The iterator can become invalid if changes are made to the repository (deletions) 
	 */ 
	class Iterator {
	public:
		/**
		* @brief Returns a pointer to the current family and advances to the next
		* @return pointer to a chromosome or NULL if we are at the end
		*/
		Pedigree *GetNext();							
		~Iterator();									///<Destruction
		void Reset();									///<Reset the iterator
		///Basic contruction
		Iterator(std::map<std::string, Pedigree*> *repository);		
	protected:	
		///<The repository we are iterating through
		std::map<std::string, Pedigree*> *repository;				
		///<The current position
		std::map<std::string, Pedigree*>::iterator position;			
		
	};

	/**
	 * @brief Write the data to the file, filename
	 */
	void Write(const char *filename);
	/**
	 * @brief Loads pedigree data from filename and returns the snp count
	 * @note eventually, we probably want to allow for different types of parsers
	 */
	uint Load(const char *filename, std::ostream& os);
	uint PostLoad(std::ostream& os);
	/**
	 * @brief builds up the records associated with each pedigree in the repository
	 */
	int BuildFamilyRecords(std::vector<FamilyRecord>& records);

	PedigreeRepository::Iterator GetIterator() { 
		return PedigreeRepository::Iterator(&pedigrees); 
	}

	/**
 	 * @brief Acquire a pointer to a pedigree. May or may not already exist in the pool
 	 */
	Pedigree *GetPedigree(const char *pedID, bool createIfNotPresent = false);

	/**
	 * @brief Acquire a pointer to an individual (may or may not exist)
	 */
	Individual *GetIndividual(const char *pedID, const char *indID, bool createIfNotPresent = false);
	
	/**
	 * @brief Describe the contents of the repository
	 */
	void GenerateReport(std::ostream &os);

	/**
	 * @brief Set up the genotypes for analysis
	 */
	bool InitializeData(GenotypeRepository& repos, Utility::Random& ran, int xvCount,  bool randomizeStatus = false);

	/**
	 * @brief devide into folds
	 */
	void InitCrossValidation(int foldCount);

	int GetPedigreeCount() { return pedigrees.size(); }

	int GetDSPCount() { return dspCount; }

	int GetLocusCount() { return snpCount; }

	void LoadDat(const char *filename);

//	static std::string statusValue;			///<What the user says is affected
	static bool verboseFoldingReport;		///<Report on the results of the xv folding

	std::vector<std::string>* GetSnpLabels() { return &snpLabels; }
protected:
	int *metaData;							///<Meta information about individuals at a given column

	//GenotypeRepository genotypeData;		///<raw genotype information to be used for analysis

	/**
	 * @brief pedigree objects, each containing the individuals that are to be analyzed
	 */
	std::map<std::string, Pedigree*>   pedigrees;


	int dspCount;							///<the number of DSPs associated with ALL pedigrees

	int snpCount;							///<The number of SNPs observed

	std::vector<std::string> snpLabels;			///<Labels for reporting (from dat file)
};


template<class gtConv>
class FamilyRepository : public PedigreeRepository {
public:
	bool ParseLine(const char *line, uint val);
protected:
	gtConv convertGenotypes;				///<Genotype Conversion	
};



template<class gtConv>
inline
bool FamilyRepository<gtConv>::ParseLine(const char *line, uint val) {
	int snpCount = 0;
	std::stringstream ss(line);
	std::string familyID, indID, dad, mom, firstOffSpringID, nextPatSibID, nextMatSibID, probandStatus;
	char status, gender;

	if (Individual::UseMerlin)
		ss>>familyID>>indID>>dad>>mom>>gender>>status;
	else
		ss>>familyID>>indID>>dad>>mom>>firstOffSpringID>>nextPatSibID>>nextMatSibID>>gender>>probandStatus>>status;

	if (familyID == "")
		return false;
	Individual *newIndividual = GetIndividual(familyID.c_str(), indID.c_str(), true);
	newIndividual->SetParentIDs(dad, mom);
	newIndividual->SetStatus(status);	
	newIndividual->SetGender(gender);
	assert(newIndividual);

	if (newIndividual->CountGenotypes() > 0) {
		std::cerr<<"An error has been encountered parsing the inputfile. Individual "<<familyID<<" x "<<indID<<" appears more than one time. Unable to continue\n";
		exit(0);
	}
	while (!ss.eof()) {
		std::string al1="", al2="";
		ss>>al1>>al2;
		
		if (al1 != "" && al2 != "") {
			newIndividual->AddGenotype(convertGenotypes(al1.c_str(), al2.c_str(), val));
			snpCount++;
		}
	}
	if (snpCount == 0)
		return true;
	if (this->snpCount > 0) {
		assert(this->snpCount == snpCount);
	}
	else
		this->snpCount = snpCount;
	return true;
}





inline
void PedigreeRepository::Iterator::Reset() {
	position = repository->begin();
}

inline
PedigreeRepository::Iterator::Iterator(
		std::map<std::string, Pedigree*> *repository) : 
				repository(repository) {
	position = repository->begin();
}



inline
PedigreeRepository::Iterator::~Iterator() { }

inline
Pedigree *PedigreeRepository::Iterator::GetNext() {
	Pedigree *current = NULL;
	if (position != repository->end()) { 
		current=position->second;
		position++;
	}
	return current;
}


}

#endif
