//
// C++ Interface: penfilemodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_STATUSMODELPENFILEMODEL_H
#define SIMULATION_STATUSMODELPENFILEMODEL_H

#include "penetrancemodel.h"



namespace Simulation {

namespace StatusModel {

/**
Used to represent penetrance models that are loaded from a file

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PenFileModel : public PenetranceModel, public Utility::AsciiParser {
public:	
	/**
	 * @brief Main Constructor
	 * @param modeID The unique ID for the model
	 * @param modelSize The number of snps in the model
	 * @param freq The frequency this model should be drawn when performing an assignment
 	 */
	PenFileModel();

	virtual ~PenFileModel();    

	/**
	 * @brief Grabs the relavent details from lines in the file
	 */
	bool ParseLine(const char *line, uint val);

	string GetModelConfiguration();

	/**
	 * @brief Load the model at the local file, filename
	 */
	virtual void Load();

	/**
	 * @brief Sets up necessary details for a given model size
	void SetModelSize(uint modelSize);
	 */

	virtual bool Init(istream& s, PoolManager *pools);

	/**
	 * @Brief Extracts filename from the stream (handling quoted filenames properly)
	 */
	//string ParseFilename(istream &s);

	/**
	 * @brief Produces basic configuration report details
	 * @param os This is the stream to which the data is written
	 * @param headerWidth This is the width of the header column. 
	 */
	void GenerateReport(std::ostream &os, uint headerWidth);

	void GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci);	

	/**
	 * @brief Set the filename to allow the program to load the data later (and not know where to load from
	 */
	void SetFilename(const char *filename);

	/**
	 * @brief Used for reporting, this just returns a user recognizable type for the model
	 */
	string Details() 		{ return string("Penetrance Model ") + filename; }

	/**	
 	 * @brief Update the allele frequencies for the model based on the contents inside the pool	
	 * This is necessary anytime the pool contents change (such as after 1 or more generational advancements) 
	 */
	virtual void Refresh(PoolManager *pools);
	int GetType() { return 2; }
protected:
	static void InvertAlleles(vector<double> &genotypes, int genotypeCount, int genotypeToInvert);
		
	/**
	 * @brief Verifies the allele falls within the allowable threshold 
	 */
	void VerifyAlleleFreq(char allele, float observedFreq);

	/**
	 * @brief Load the penetrance table and other details from the file	
	 */
	virtual void Load(const char *filename);


	uint modelSize;						///<Local size of the disease model
	std::string filename;				///<Where the model details are to be read/written
	float freqThreshold;				///<Allowable variation of allele frequencies
	map<char, float> alleleFreq;		///<Lookup for allele frequencies
};

}

}

#endif
