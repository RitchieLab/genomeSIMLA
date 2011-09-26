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
#ifndef SIMULATION_STATUSMODELCONTINOUSMODEL_H
#define SIMULATION_STATUSMODELCONTINOUSMODEL_H

#include "penetrancemodel.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/shared_ptr.hpp>
#include "utility/random.h"
#include "diseasemodel.h"
#include "utility/lineparser.h"
#include "simulation/poolmanager.h"


namespace Simulation {

namespace StatusModel {


#define LABEL_GRANDMEAN "GRAND_MEAN"
#define LABEL_LOWER_THRESHOLD "MIN_THRESHOLD"
#define LABEL_UPPER_THRESHOLD "MAX_THRESHOLD"
/**
 * @brief Simple struct to aggregate the distribution with the mean/std-dev
 */
struct NormalDistribution {
	/**
	 * @brief Basic constructor
	 */
	NormalDistribution(int genotype = -1, float mean = 0.0, float sigma=0.0) : 
			genotype(genotype), distribution(mean, sigma) { }


	void InitDistribution(float mean, float sigma) {
		distribution = boost::normal_distribution<float>(mean, sigma);
	}
	/**
	 * @brief Perform the drawn (pass it a copy of the appropriate generator)
	 */
	float GetOutcome(Random& rnd) {

		return distribution(rnd);
	}

	/** Just book keeping-ID one's self */
	int genotype;							

	/** Used to draw the random numbers from the appropriate stream */
	boost::normal_distribution<float> distribution;
};

typedef boost::shared_ptr<NormalDistribution> DistributionPointer;

/**
Used to represent penetrance models that are loaded from a file

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ContinuousModel : public PenetranceModel, public Utility::AsciiParser {
public:	
	/**
	 * @brief Main Constructor
	 * @param modeID The unique ID for the model
	 * @param modelSize The number of snps in the model
	 * @param freq The frequency this model should be drawn when performing an assignment
	 */
	ContinuousModel();

	virtual ~ContinuousModel();    

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
	string Details() 		{ return string("Continous Outcome Model ") + filename; }

	/**	
 	 * @brief Update the allele frequencies for the model based on the contents inside the pool	
	 * This is necessary anytime the pool contents change (such as after 1 or more generational advancements) 
	 */
	virtual void Refresh(PoolManager *pools);

	/**
	 * Performs a quick evaluation of the individual's genetic makeup and determines status randomly
	 * @return 012 where 0 is lower tail, 1 is upper tail and 2 is interior 
	 * 
	 */
	virtual int GetStatus(std::vector<uint>& genotypes, float& outcome);

	void ReportTable(ostream& os);

	void AddCell(uint idx, double mean, double stddev);

	/**
	 * @brief Used to split the distribution into tails in order to allow the dataset to be proportionally selected for
	 */
	void StdDeviationsFromMean(float& count);

	/**
	 * @brief returns the group associated with the status, -1 if inviable
 	 */
	int DrawPopulationStatus(float& status);
	int DetermineStatus(float outcome);
	/**
	 * @brief Indicates that the model is continuous (or discrete)
	 */
	virtual bool IsContinuous() { return true;}

	void GetTailBounds(float& lower, float& upper) { lower=lowerTailMax; upper=upperTailMin; }
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


	uint modelSize;							///<Local size of the disease model
	std::string filename;					///<Where the model details are to be read/written
	float freqThreshold;					///<Allowable variation of allele frequencies
	map<char, float> alleleFreq;			///<Lookup for allele frequencies
	vector<DistributionPointer> dist;		///<This is the array for converting genotypes to distributions
	NormalDistribution globalDistribution;	///<The population mean/std-dev
	float lowerTailMax;						///<max bound for the lower tail
	float upperTailMin;						///<Min bound for the upper tail
	float numDev;							///<The number of deviations from the mean to separate the tails
	int totalQueries;						///<Used to calculate the pevalence 
	
	//An individual whose outcome is 'unviable' should be discarded
	int modelThreshold;						///<0 no threshold, 1 lower, 2 upper, 3 both
	float minThreshold;						///<Lower bound for viable outcomes
	float maxThreshold;						///<Upper bound for viable outcomes
};

}

}

#endif
