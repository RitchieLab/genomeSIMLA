//
// C++ Interface: pedigree
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTPEDIGREE_H
#define MDRPDTPEDIGREE_H

#include <map>
#include <string>
#include <iostream>
#include "individual.h"
#include "sibship.h"

namespace MdrPDT {



/**
 * @brief Organization structure used to properly reference each of the members of the family
 * @note It is imperative that all individuals be created through the pedigree they are to be
 * 		assigned. This insures that they are properly assigned to their sibships and all.

	@author Eric Torstenson
*/
class Pedigree{
public:
    Pedigree(const char *id, int pedID);
    ~Pedigree();
	/**
	 * @brief This allows basic iteration over the members of the pedigree. 
	 * The iterator can become invalid if changes are made to the repository (deletions) 
	 */ 
	class Iterator {
	public:
		/**
		* @brief Returns a pointer to the current individual and advances to the next
		* @return pointer to a individual or NULL if we are at the end
		*/
		Individual *GetNext();							
		~Iterator();									///<Destruction
		void Reset();									///<Reset the iterator
		///Basic contruction
		Iterator(std::map<std::string, Individual *> *repository);		
	protected:	
		///<The repository we are iterating through
		std::map<std::string, Individual *> *repository;				
		///<The current position
		std::map<std::string, Individual *>::iterator position;			
		
	};
	
	/**	
 	 * @brief Returns a pointer to an individual, or NULL
	 * @note Depending on createIfNotPResent, this function might allocate new memory
	 */
	Individual *GetMember(const char *id, bool createIfNotPresent = false);

	Pedigree::Iterator GetIterator() { 
		return Pedigree::Iterator(&members); 
	}

	/** 
	 * Some things require some postwork to be performed. 
	 * Some details to be performed include:
	 * 		* Define the pedigreeID...an integer value associated with valid pedigrees
	 * 		* Create the sibships
	 * 		* Expand "virtual" siblings
	 */
	void PostLoad(std::ostream& os);						

	/**
	 * @brief Present some details found in the dataset, and possibly some basic settings
	 */
	void GenerateReport(std::ostream &os);
	
	/**
	 * @brief Return the number of DSPs associated with the pedigree
	 */
	int GetDSPCount();

	/**
	 * @brief Returns the id as it was seen in the dataset
	 */
	std::string ID() { return id; }

	/**
	 * @brief Writes genotype information to the grid for all DSPs
	 */
	int WriteGenotypeData(char *data, int *meta, char *folds, int offset, Utility::Random& gen, int xvSlice, int stride, bool permute = false);

	void Write(std::ostream& os);

	bool DropFromAnalysis() { return dropFromAnalysis; }
	void DropFromAnalysis(bool doDrop) { dropFromAnalysis = doDrop; }

	/**
	 * @brief Fills in the rest of the familial details once all individuals have been loaded
	 */
	void Reconcile();
	/**
	 * @brief Returns the number of members are part of the pedigree
	 */
	int GetMemberCount();
protected:
	/**
	 * @brief Purge the members of the various storage devices
	 */
	void PurgeMembers();

	/**
	 * @brief Initialize the sibship structure. This will create the existing sibships based on the dataset
	 */
	void InitSibships();
	
	/**
	 * @brief Expand the sibships according the trio/expand all. 
	 * @note This should be done ONLY one time, and immediately after the data has been loaded
	 */
	void ExpandSibships(std::ostream& os);

	std::string id;
	std::map<std::string, Individual *> members;

	/**
	 * @brief unique number which will be used in the construction of the D-Statistic, which reside inside an array
	 */
	int pedigreeID;

	//This map doesn't hold anything unique as far as memory goes....it doens't need to be "purged", though, it should be cleared when the members are cleared
	std::map<std::string, Sibship > sibships;

	/**
	 * @brief cached DSP count. returns -1 if there are no sibships associated with the family
	 */
	int dspCount;

	bool dropFromAnalysis;

};


inline
void Pedigree::Iterator::Reset() {
	position = repository->begin();
}

inline
Pedigree::Iterator::Iterator(
		std::map<std::string, Individual*> *repository) : 
				repository(repository) {
	position = repository->begin();
}



inline
Pedigree::Iterator::~Iterator() { }

inline
Individual *Pedigree::Iterator::GetNext() {
	Individual *current = NULL;
	if (position != repository->end()) { 
		current=position->second;
		position++;
	}
	return current;
}


}

#endif
