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


// blockpool.h

/**
 * @brief Maintains the details associated with a single block of dna data
 */

#ifndef __CHROMPOOL_H__
#define __CHROMPOOL_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "chromosome.h"
#include "individual.h"
#include "growthrate.h"
#include "distmappingfn.h"
//#include "population.h"

namespace Simulation {


/**
 * @brief Manages the various expressions of a pool of chromosomal data
 * ChromPool objects can be intialized and advanced through generations using
 * "sexual reproduction" in which two expressions are plucked from the pool (with replacement) 
 * and crossed. The various points of cross over are determined by drawing against the locus
 * recombination point. 
 * Blocks are produced with different recombination fraction ranges from the overall chromosomal
 * values - This allows the user to control how likely the LD between blocks is. 
 * ChromPool Objects are totally independant of one another and ultimately all operations on them
 * will be threadsafe.
 */


class ChromPool {

public:
	struct BlockDefinition {
		uint minSnpCount;				///<Minimum number of snps
		uint maxSnpCount;				///<Maximum number of snps
		float minBlckRecomb;			///<Recombination factor for the second locus
		float maxBlckRecomb;			///<Recombination factor for first locus
		float minSnpRecomb;				///<min Recombination factor for the rest of the snps
		float maxSnpRecomb;				///<Max recombination fraction for the rest of the SNPS
		//float strength;					///<Strength of the block's front end
		float frequency;				///<Likelihood of being drawn
		uint id;						///<Used to report which block we are talking about

		BlockDefinition() : minSnpCount(0), maxSnpCount(100), minBlckRecomb(0.01), maxBlckRecomb(0.05), minSnpRecomb(0.00001), maxSnpRecomb(0.00005), frequency(0.0), id(0) { }

		BlockDefinition(uint min, uint max, float minBlckRecomb, float maxBlckRecomb, float minSnpRecomb, float maxSnpRecomb, float freq, uint id) : 
			minSnpCount(min), maxSnpCount(max), 
			minBlckRecomb(minBlckRecomb), maxBlckRecomb(maxBlckRecomb),
			minSnpRecomb(minSnpRecomb), maxSnpRecomb(maxSnpRecomb), frequency(freq), id(id) {}

		/**
		 * @brief let the definition determine these values
		 */
		uint GetSnpCount();
		float GetBlockRecombination();
		float GetSnpRecombination();
	} defaultBlock;	
	/**
	 * @brief basic construction
	 * @param chromID This is used to identify which chromosome this block is assigned to
	 * @param blockID Used to uniquely id the local block
	 */
    ChromPool(uint chromID, uint blockCount, float minRecomb, float maxRecomb);
	ChromPool(uint chromID, const char *locFilename);
	~ChromPool() {}

	/**
	 * Save in phased format ready for use in haploview
	 */
	void SaveAsPhased(std::ostream &os);	

	/**
	 * Save in phased format starting with first, and ending in last
	 */
	uint SaveAsPhased(std::ostream& os, uint start, uint count);
	/**
	 * @brief updates allele frequencies.
	 * @note Call this after any changes have been made to the pool (this is totally unrelated to data sets)
	 */
	void CalculateAlleleFrequencies();

	/**
	 * @brief Load loci from a file or input stream. This should be done before calling LoadPhased
	 */
	uint LoadLoci(istream& data);
	

	/**
	 * @brief Save loci to stream
	 */
	void SaveLoci(ostream& os);

	/**
	 * @brief Returns true if there is enough information to produce map info files for haploview charts
	 */
	bool UseMapInfoFile();

	/**
	 * @brief Loads the phased data into the local pool. 
	 * @param os The stream containing the data
	 * @param currGen The number of generations already ran
	 */
	void LoadPhased(std::istream &os, uint currGen);

	/**
	 * @brief Add a potential block of loci the chromosome
	 * @param minSnpCount The minimum number of snps 
	 * @param maxSnpCount The the maximum number of snps
	 * @param minRec The minimum recombination Fraction between one locus and and the next one
	 * @param maxRec The max recombination fraction between one locus and the next one
	 * @param frequency The frequency with which this block should be drawn
	 */
	uint DefineBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax, float frequency);
		
	/**
	 * @brief Defines the default block
	 * @param minSnpCount The minimum number of snps 
	 * @param maxSnpCount The the maximum number of snps
	 * @param minRec The minimum recombination Fraction between one locus and and the next one
	 * @param maxRec The max recombination fraction between one locus and the next one
	 * @note Default blocks are used (1 per chromosome) to fill in when the user's block definition frequencies don't reach 1.0
	 */
	void DefineDefaultBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax);

	void DefineDefaultBlock(ChromPool::BlockDefinition& block);

	/**
	 * @brief Initialize the loci associated with the local chromosome
	 * @note Call this after defining all blocks and before population generation. LoadLoci takes care of this when loading from a file
	 */
	uint InitLoci();

	/**
	 * @brief returns the recombination fraction for a given locus
	 */
	double GetRecombinationFraction(uint locus);


	/**
	 * Builds up the pool of blocks
	 */
	void BuildPool(uint expressionCount);

	/**
	 * @brief Forces the allele frequency to be al1 and al2 at locus, locusID
	 * @param locusID The index (0 based) of the locus of interest
	 * @param al1 The first allele's frequency
	 * @param al2 The frequency of the second allele
	 */
	bool ForceAlleleFrequency(uint locusID, float al1, float al2);

	/**
	 * @brief Draws a unique individual from the pool, guaranteeing that the individual is unique
	 * @note When Returning an individual to the pool, they can be drawn later. 
	 */
	void DrawIndividual(Individual *);

	/**
	 * @brief returns the error rate for a given locus
	 */
	float GetErrorRate(uint locusID);

	/**
	 * @brief Returns the id associated with the local chromosome
	 */
	uint GetID() { return chromID; }

	/**
	 * @brief Returns the number of loci associated with the local chromosome
	 */
	uint GetLociCount() { return lociCount; }

	/**
	 * @brief Returns the number of individual expressions of a given chromosome inside the pool
	 */
	uint GetExpressionCount() { return pool.size(); }
	/**
	 * Performs generational mating/replacement
	 * @param genCount The number of generations to move forward. 
	 * @param endSize The desired size after all generations have been completed. 0 means no change
	 * @return The current generation after we are finished
	 */
	uint AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f);
	

	/**
	 * @brief Empties the pool and resets the generation count
	 */
    void Clear();

	/**
	 * @brief Allows an individual to be redrawn
	 * @note Only use this when an individual is known to not exist in any of the samples. 
	 */
	void ReturnIndividual(Individual* ind);

	static float defFre1;				///<Generic value used in the initial generation of loci
	static float defFre2;				///<Generic value used in the initial generation of loci
	static float errRate;				//
	static bool randomizeAlleleFreq;	///<Set by configuration to determine how to set up the loci
	
	//Growth rate details
	static float minGrowthRate;			///<Min percentage of current population (can be +/-)
	static float maxGrowthRate;			///<Max percentage of current population (can be +/-)
	static bool allowPopulationGrowth;	///<Turn growth on or off

	static string javaPath;
	static string haploviewPath;
	static string haploviewSettings;	

	static DistMappingFn *mappingFn;	///<The mapping function to be used to convert recomb fraction to bases

	typedef std::vector<Chromosome> PoolType;
	typedef std::map<std::string, uint> IndLookupType;


	/**
	 * @brief Calculates the allele frequencies
	 */
	bool GetAlleleFrequencies(uint locus, float &af1, float &af2);

	/**
	 * @brief Convenience overload for complete chromosome information
 	 * @param os The stream to be written to
	 * @param f The function used to determine the distance
	 * @note This just iterates from locus 0 to locus lociCount
	 */
	void WriteMarkerInfo(std::ostream& os, DistMappingFn *f=NULL);

	/**
	 * @brief Write the marker information file to a stream
	 * @param os The stream to be written to
	 * @param first The first locus to be written
	 * @param last The last locus to be written
	 * @param f The Mapping function to be used. If NULL, we will use the local fn
	 * @note Basically, work through the loci from first to last. If it is important
	 * f is NULL, the local function is used and the position will be reset to 0. 
	 */
	void WriteMarkerInfo(std::ostream& os, uint first, uint last, DistMappingFn *f=NULL);

private:
	bool frequenciesDetermined;
	///Used to determine the number of columns are in a line
	uint CountColumns(const char *line);				
	

	///<Simple cheat to check if we've created this individual before
	IndLookupType usedIndividuals;	

    PoolType pool;						///<The actual pool
	uint chromID;						///<Indicate which chromosome the block is associated with
  	uint generationCount;				///<How many generations lead us to this point
	uint lociCount;						///<Number of loci associated with the local pool
	uint expressionCount;				///<Indicate the size of the pool

	float minRecomb;					///<Minimum Locus recombination fraction
	float maxRecomb;					///<Max locus recombination fraction
	uint blockCount;					///<Used to initialize the chromosome's loci

	LocusArray loci;					///<Description of the various loci associated with the chromosomes
	vector<BlockDefinition> blockPrototypes;
	BlockDefinition &DrawBlockDefinition();

	string locSource;					///<Used for a chromosome that is to be based on a source file
	bool positionsDefined;				///<Record that we have real positions defined
};


inline
uint ChromPool::BlockDefinition::GetSnpCount() {
	//Figure out how many loci there are locally
	int diff = maxSnpCount - minSnpCount;
	return minSnpCount + Utility::Random::globalGenerator(diff);
}

inline
float ChromPool::BlockDefinition::GetBlockRecombination() {
	//First one we use the block strength variation from 0.5
	float rDiff = maxBlckRecomb - minBlckRecomb;
	return minBlckRecomb + (Utility::Random::globalGenerator.drand() * rDiff);
}

inline
float ChromPool::BlockDefinition::GetSnpRecombination() {
	float rDiff = maxSnpRecomb - minSnpRecomb;
	return minSnpRecomb + (Utility::Random::globalGenerator.drand() * rDiff);
}

}

#endif
