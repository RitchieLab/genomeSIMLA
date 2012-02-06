//
// C++ Interface: diseasemodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONDISEASEMODEL_H
#define SIMULATIONDISEASEMODEL_H


//#define DEBUG_PRODUCTION 1

#include "utility/types.h"
#include "diseaselocus.h"

namespace Simulation {
namespace StatusModel {

using namespace Utility;

/**
 * @brief Base class for disease models. 
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class DiseaseModel  {
public:
	/**
	 * @brief construction
	 * @param modelID The unique ID for the model	
 	 * @param prob The probability an individual will be affected with this type
	 */
	DiseaseModel(uint modelID = 0) : modelID(modelID) {}
	virtual ~DiseaseModel() {}

	/**
	 * @brief Returns probability associated with the local model
	 */
	virtual float GetProbability() { return probability; } 

	/**
	 * @brief Returns the unique ID for the model 
	 */
	virtual uint GetModelID(uint modelID) { return modelID; }

	/**
	 * Performs a quick evaluation of the individual's genetic makeup and determines status randomly
	 * 
	 */
	virtual int GetStatus(std::vector<uint>& genotypes, float& outcome) = 0;

	/**
	 * @brief Returns the number of loci associated with the model	
	 */
	virtual size_t GetModelSize() = 0;


	/**
	 * @brief Report the details of this model to the stream
	 * @param os Stream to be written to
	 * @param headerWidth The number of spaces on each line before we print anything
	 */
	virtual void GenerateReport(std::ostream &os, uint headerWidth) = 0;
	/**
	 * @brief Each model will describe itself for reporting purposes
	 */
	virtual std::string Details() = 0;

	/**
	 * @brief Add a single locus and penetrance to the table
	 * @param chromID The chromosome on which the snp is found
	 * @param locus The snp locus being described
	 * @param minorAlleleFreq  
	 * @return Returns the size of the model thus far
	 */
	void AddDiseaseLoci(const char *label, uint chromID, uint locus, float minorAlleleFreq);	
	
	/**
	 * @Brief Returns the locus at a given position (pos must be < GetModelSize() )
	 */
	DiseaseLocus GetLocus(uint pos);

	/**
	 * @brief Returns a bitset containing true for each chromosome Index that is linked to the model
	 * @param chromCount The number of chromosomes in the whole genome
	 */
	BitSetType AssociatedChromosomes(uint chromCount);

	/**
	 * @brief Returns the various loci found on a given chromosome
	 * @param chromID the chromosome index to search for
	 * @param relLoc The vector into which the indices will be written
	 * @return t/f indicating whether or not anything was found
	 */
	bool GetDiseaseLoci(int chromID, vector<uint> &relLoc);

	/**
	 * @brief Returns vector containing each of the loci associated with the model
	 */
	ModelLociArray GetDiseaseLoci() { return loci; }

	virtual int GetType() { return 0; }
protected:
	uint modelID;						///<Unique ID for a given model
	float probability;					///<This is used to determine percentage of a population to be assed by this model
	ModelLociArray loci;				///<Used to store locus IDs associated with model


};

}
}

#endif
