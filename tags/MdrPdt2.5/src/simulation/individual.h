//
// C++ Interface: individual
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONINDIVIDUAL_H
#define SIMULATIONINDIVIDUAL_H

#include <string>
#include <vector>
#include <map>
#include "chromosome.h"
#include "penetrancemodel.h"
#include "utility/utility.h"

namespace Simulation {

using namespace std;

/**
base class for a person

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Individual{
public:
	/**
	 * @brief Construction for arrays: Usage requires a call to Init() afterward to get everything set up properly
	 */
	Individual();
    Individual(uint id, uint pedID, uint chrCount=1);
    virtual ~Individual();

	/**
	 * @brief Initialize the basic components (if using default constructor)
	 */
	void Init(uint id, uint pedID, uint chrCount=1);

	uint GetID() { return id; }

	uint GetPedigreeID() { return pedID; }

	/**
	 * Return the pool ID associated with the maternal and paternal Chromosomes. 
	 * This value is used to help avoid duplicate individuals being drawn.
	 */
	std::string GetPoolID(uint chromeID);

	/**
	 * Set an ID associated with a given chromosome
	 */
	void SetPoolID(uint chromeID, const char *id);

	/**	
	 * Set the chromosomal data
	 */
	void SetChromosomalData(uint chromeID, Chromosome &p, Chromosome &m);
	void SetChromosomalData(Chromosome *p, Chromosome *m);
	/**
	 * @brief Pass a model to the individual to allow it to determine the status
	 * @param model This is the model that is to be used in assigning status
	 * @return T/F indicating affected status
	 */
	virtual bool ApplyStatus(PenetranceModel* model);

	/**
	 * @brief Returns the genotype for a given locus
	 */	
	virtual int GetGenotype(uint chromID, uint locus);

	virtual int GetSnpCount();

	/**
	 * @brief flip one of the alleles according the loci direction
	 */
	virtual bool ChangeGenotype(uint chrID, uint locus, int errorDir);

	virtual bool ClearLocus(uint chrID, uint loc);
	/**
	 * @brief Setup the affected value
	 */
	virtual void SetStatus(bool isAffected);
	
	/**
	 * @brief Used to determine if an individual is affected
	 */
	virtual bool IsAffected();

	/**
	 * @brief this only matters if you want valid pedigrees to be written
	 */
	void SetPedigreeMeta(uint patID, uint MatID);

	void WritePedigree(ostream& os, uint *genotypeCounts);
	void WritePhased(ostream& os, uint indID);
	void WriteMDR(ostream& os, uint *genotypeCounts);

	/**
	 * @brief Produce a brand new individual with genetic information from both the mother and the father (the local is considered mommy)
	 */
	Individual *Cross(Individual *other, uint id);
	
	/**
	 * @brief Produce a new chromosome array based on a crossing of the local chromosomal data
	 */
	Chromosome *Cross();

	bool DoIncludeInDataset();
	void DoIncludeInDataset(bool doInclude);

	/**
	 * @brief Ordinarily, we want our data with 10 header columns. However, some apps require 6. False sets it to 6.
	 */
	static bool StandardPedigreeHeader;
protected:

	/**
	 * Pedigree information
	 */
	uint id;								///<The individual id
	uint pedID;								///<pedigree id
	uint patID;								///<Paternal id
	uint matID;								///<Maternal id
	bool includeInDataset;					///<Allow the system to mark individuals to not be written to datasets

    unsigned int status;					///<Overall status to be reported in the file
	Chromosome *matChrom;					///<Mother's contribution
	Chromosome *patChrom;					///<Father's contribution
	uint chromosomeCount;					///<the number of chromosomes associated with the individual
	string *genePoolIDs;					///<Used so we can return an individual to the pool to be redrawn later if it couldn't be used for some reason

	Utility::BitSetType *missingData;



};


}

#endif
