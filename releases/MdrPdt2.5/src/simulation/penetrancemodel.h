//Model.h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////// 
//
// Class contains penetrance tables and 
// list of loci associated with the model.
// The list of loci is kept as the locus indexes
// as kept in the Chromosome static data.
//
/////////////////////////////////////////////////////////////////////


#ifndef __BLANKMODEL_H__
#define __BLANKMODEL_H__

#include <string>
#include <iostream>
#include <vector>
#include "utility/utility.h"
#include "diseasemodel.h"


/*******************smod file***********************/
#define LABEL_DISEASELOCI "DISEASELOCI"
#define LABEL_PENTABLE "PENTABLE"


namespace Simulation {



class PenetranceModel : public DiseaseModel, public Utility::AsciiParser {
public:	
	/**
	 * @brief Main Constructor
	 * @param modeID The unique ID for the model
	 * @param modelSize The number of snps in the model
	 * @param freq The frequency this model should be drawn when performing an assignment
 	 */
	PenetranceModel(uint modelID, uint modelSize, float prob);

	/**
	 * @brief Constructor for use when loading from a file
	 * @param modeID The unique ID for the model
	 * @param freq Frequency this model should be drwan when performing an assignment
	 */
	PenetranceModel(uint modelID, float prob);

	virtual ~PenetranceModel();    
	
	/**
	 * @brief Add a single locus and penetrance to the table
	 * @param chromID The chromosome on which the snp is found
	 * @param locus The snp locus being described
	 * @param penetrance 
	 * @return Returns the size of the model thus far
	 */
	void AddDiseaseLoci(uint idx, uint chromID, uint locus);

	/**	
 	 * @brief We are assuming that idx will be 0 based and ordered correctly
	 * AABBCCDD, AABBCCDd, AABBCCdd, AABBCcDD, AABBCcDd, etc....
	 */
	void AddPenetrance(uint idx, double penetrance);
		
	/**
	 * @brief Returns the locus number at pos
	 */
	uint GetLocus(uint pos);

	/**
	 * @brief Returns the chromosome ID at pos
	 */
	uint GetChromID(uint pos);


	/**
	 * @brief Returns the number of loci associated with the model
	 */
	uint GetModelSize();

	/**
	 * @brief Grabs the relavent details from lines in the file
	 */
	bool ParseLine(const char *line, uint val);

	/**
	 * @brief Converts the letters, AABb, etc... into a numerical index
	 */
	uint GetGenotypeIdx(const char *genotype, uint lociCount);

	/**
	 * Performs a quick evaluation of the individual's genetic makeup and determines status randomly
	 * 
	 */
	virtual bool IsAffected(std::vector<uint>& genotypes);

	/**
	 * @brief Load the penetrance table and other details from the file	
	 */
	virtual void Load(const char *filename);

	/**
	 * @brief Load the model at the local file, filename
	 */
	virtual void Load();

	/**
	 * @brief Produces basic configuration report details
	 * @param os This is the stream to which the data is written
	 * @param headerWidth This is the width of the header column. 
	 */
	void GenerateReport(std::ostream &os, uint headerWidth);
	
	/**
	 * @brief Set the filename to allow the program to load the data later (and not know where to load from
	 */
	void SetFilename(const char *filename);

private:


/**
 * @brief Help keep up with chromosome and locus indices
 */
struct LocusType {
	uint chromosomeID;					
	uint locusID;	

	LocusType(uint c, uint l) : chromosomeID(c), locusID(l) {} 
	LocusType() : chromosomeID(0), locusID(0) {}
};
	/**
	 * @brief Configure the disease loci based on input from the model file
	 */
	void SetupDiseaseLoci(const char *line);

	/**
	 * @brief Used in determining the position with the model's penetrance table based on genotype
	 */
	uint GetMultiplier(uint genotype, uint position);


	uint penCount;						///<Size of the penList array
    double *penList;					///<The array of penetrances for each of the loci
	uint maxPenCount;					///<The maximum number of penetrances possible for a given sized model	
 
   	LocusType *modelLoci;				///<Array of loci associated with the local model
	uint modelSize;						///<Local size of the disease model

	bool isConfigured;					///<Determines if the model is properly configured
	std::string filename;				///<Where the model details are to be read/written
};



}


#endif 
