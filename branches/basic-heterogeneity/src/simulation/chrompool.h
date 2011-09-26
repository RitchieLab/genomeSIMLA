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

#include "minchrom.h"
#include "chromosome.h"
#include "growthrate.h"
#include "distmappingfn.h"
#include "individual.h"
#include "blocklistnodefourgammetes.h"
#include <map>
#include "utility/rbtree.h"
#include "locreport.h"
#include "utility/exception.h"
#include <pthread.h>
//#include "haplotypeblock.h"
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

typedef std::map<std::string, Locus> LocusMap;

#define LOCKPOOL 	pthread_mutex_lock(&poolLock)
#define UNLOCKPOOL 	pthread_mutex_unlock(&poolLock)

#define LOCKREPORT  pthread_mutex_lock(&reportLock)
#define UNLOCKREPORT pthread_mutex_unlock(&reportLock)

class ChromPool {
public:

	typedef enum ChromosomeType { BlockBased, FileBased } ChromosomeType;

	struct BlockDefinition {
		uint minSnpCount;				///<Minimum number of snps
		uint maxSnpCount;				///<Maximum number of snps
		float minBlckMap;			///<Recombination factor for the second locus
		float maxBlckMap;			///<Recombination factor for first locus
		float minSnpMap;				///<min Recombination factor for the rest of the snps
		float maxSnpMap;				///<Max recombination fraction for the rest of the SNPS
		//float strength;					///<Strength of the block's front end
		float frequency;				///<Likelihood of being drawn
		uint id;						///<Used to report which block we are talking about

		BlockDefinition() : minSnpCount(5), maxSnpCount(100), minBlckMap((float)0.0001), maxBlckMap((float)0.0005), minSnpMap((float)0.00001), maxSnpMap((float)0.00005), frequency((float)0.0), id(0) { }

		BlockDefinition(uint min, uint max, float minBlckMap, float maxBlckMap, float minSnpMap, float maxSnpMap, float freq, uint id) : 
			minSnpCount(min), maxSnpCount(max), 
			minBlckMap(minBlckMap), maxBlckMap(maxBlckMap),
			minSnpMap(minSnpMap), maxSnpMap(maxSnpMap), frequency(freq), id(id) {}


		/**
		 * @brief Create random values (for GUI)
		 */
		void Randomize();

		/**
		 * @brief let the definition determine these values
		 */
		uint GetSnpCount();			

		/**
		 * @brief Returns a random recombination value based on the min/max block recombination
	 	 */
		float GetBlockMapDistance();

		/**
		 * @brief Returns a random recombination value based on the min/max snp recombination
	 	 */
		float GetSnpMapDistance();


		void WriteConfiguration(ostream &str);
	} defaultBlock;	
	/**
	 * @brief basic construction
	 * @param chromID This is used to identify which chromosome this block is assigned to
	 * @param blockID Used to uniquely id the local block
	 */
    ChromPool(uint chromID, uint blockCount);
	ChromPool(uint chromID, const char *locFilename);
	ChromPool(uint chromID, const char *locfilename, const char *gridFilename);
	ChromPool(uint chromID, const char *seedFilename, const char *locFilename, const char *fmt, uint headerCount);
	~ChromPool() { if (summaryReport) { delete summaryReport; } }

	void SetBlockCount(uint bc) {	blockCount = bc; }

	void ValidateBlocks();

	/**
	 * @brief Save in phased format ready for use in haploview
	 * @param os Stream to write to
	 * @param minThreshold The minimum allele frequency to be written (use zero if you want all written)
	 */
	void SaveAsPhased(std::ostream &os, float minThreshold);	

	/**
	 * Save in phased format starting with first, and ending in last
	 * @param os Stream to write to
	 * @param minThreshold The minimum allele frequency to be written (use zero if you want all written)
	 * @param start The first loci to be written
	 * @param count the number of loci to be written to the stream 
	 */
	uint SaveAsPhased(std::ostream& os, float minThreshold, uint start, uint count);

	void SamplePhased(std::ostream& os, uint firstExpression, uint expCount, uint firstSnp, uint snpCount);
	
	/**
	 * @brief updates allele frequencies.
	 * @param maxindividualCount Default value of 0 means all "individuals"
	 * @note Call this after any changes have been made to the pool (this is totally unrelated to data sets)
	 */
	size_t CalculateAlleleFrequencies(uint maxIndividualCount = 0);

	uint LoadLoci(const char *project, LocusMap& lmap, bool forceLoad = true);

	/**
	 * @brief Load loci from a file or input stream. This should be done before calling LoadPhased
	 * @param data The stream where the data can be found
	 * @param forceLoad Will dump existing loci in leiu of the contents of the file. 
	 * @note If this is just disk caching...and the loci are already OK, set forceLoad to false
	 */
	uint LoadLoci(istream& data, LocusMap& lmap, bool forceLoad);
	

	string GenerateLocusFilename(const char *project);

	void SaveLoci (const char *project, float minThreshold = 0.0);

	/**
	 * @brief Save loci to file
	 */
//	void SaveLoci(const char *filename);
	/**
	 * @brief Save loci to stream
	 */
	void SaveLoci(ostream& os, float minThreshold = 0.0);

	/**
	 * @brief Returns true if there is enough information to produce map info files for haploview charts
	 */
	bool UseMapInfoFile();

	/**
	 * @brief Loads the phased data into the local pool. 
	 * @param os The stream containing the data
	 * @param currGen The number of generations already ran
	 * @param headerCount The number of headers (if it were from standard pedigree datasets)
	 * @param forceLoad If true, the pool will be dumped and replaced with the contents of the file. 
	 * If forceLoad=false and the pool is already populated, we won't read in the contents of the file.
	 */
	void LoadPhased(std::istream &os, uint currGen, uint headerCount, bool forceLoad);// = false);

	/**
	 * @brief Add a potential block of loci the chromosome
	 * @param minSnpCount The minimum number of snps 
	 * @param maxSnpCount The the maximum number of snps
	 * @param minRec The minimum map distance between one locus and and the next one
	 * @param maxRec The max map distance between one locus and the next one
	 * @param frequency The frequency with which this block should be drawn
	 */
	uint DefineBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax, float frequency);
		
	/**
	 * @brief Defines the default block
	 * @param minSnpCount The minimum number of snps 
	 * @param maxSnpCount The the maximum number of snps
	 * @param minRec The map distance Fraction between one locus and and the next one
	 * @param maxRec The map distance fraction between one locus and the next one
	 * @note Default blocks are used (1 per chromosome) to fill in when the user's block definition frequencies don't reach 1.0
	 */
	void DefineDefaultBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax);

	void DefineDefaultBlock(ChromPool::BlockDefinition& block);

	/**
	 * @brief Initialize the loci associated with the local chromosome
	 * @note Call this after defining all blocks and before population generation. LoadLoci takes care of this when loading from a file
	 */
	uint InitLoci(uint starterID, LocusMap &locusMap);

	/**
	 * @brief returns the recombination fraction for a given locus
	 */
	double GetMapDistance(uint locus);


	/**
	 * @brief Commits the name of the stylesheet tot he stream. 
	 * @TODO eventually, we will want to write the default style sheet if the indicated file doesn't exist
	 */
	void WriteStyleSheetDetails(ostream &os);


	/**
	 * @brief Quickly print n genotypes from the first m expressions for debugging purposes
	 */
	void SamplePool(size_t expCount, size_t genotypeCount);

	/**
	 * @brief Writes the LD information to file based on the user's preferences
	 * @param filename This is the base of the filename to be used.
	 * @param thresh the threshold to be used for excluding loci based on minor allele frequency
	 */
	void WriteLdData(const char *filename, ostream &summary, Visualization::LocusReport &locReport, float thresh = 0.0001);

	/**
	 * Builds up the pool of blocks
	 */
	void BuildPool(uint expressionCount);
/*	void BuildPool_Eden(uint expressionCount);
	void BuildPool_AdamEve(uint expressionCount);
*/	/**
	 * @brief Forces the allele frequency to be al1 and al2 at locus, locusID
	 * @param locusID The index (0 based) of the locus of interest
	 * @param al1 The first allele's frequency
	 * @param al2 The frequency of the second allele
	 */
	bool ForceAlleleFrequency(uint locusID, float al1, float al2, LocusMap &locusMap);

	/**
	 * @brief Draws a unique individual from the pool, guaranteeing that the individual is unique
	 * @note When Returning an individual to the pool, they can be drawn later. 
	 */
	void DrawIndividual(Individual *);
	void ApplyChromosomeIDs(Individual *ind);
	void ApplyPresentGenotypes(Individual *ind, bool modelLociOnly);
	/**
	 * @brief returns the error rate for a given locus
	 */
//	float GetErrorRate(uint locusID);

	/**
	 * @brief Returns the id associated with the local chromosome
	 */
	uint GetID() { return chromID; }
	
	void SetID(uint id) { chromID=id; }

	string GetLabel() { return label; }
	void SetLabel(const char *lbl) { label=lbl; }

	void SetGeneration(uint generation) { generationCount = generation; reloadLocusFile=true;}

	/**
	 * @brief Returns the number of loci associated with the local chromosome
	 */
	uint GetLociCount() { return lociCount; }

	LocusArray &GetLoci() { return loci; }

	string GetLocSource() { return locSource; }

	void SetLocSource(const char *src) { locSource = src; }

	/**
	 * @brief Returns the number of individual expressions of a given chromosome inside the pool
	 */
	size_t GetExpressionCount() { 		
		size_t size = pool.size(); 
		if (size == 0) 
			size = minPool.size(); 
		return size;
	}
	/**
	 * Performs generational mating/replacement
	 * @param genCount The number of generations to move forward. 
	 * @param endSize The desired size after all generations have been completed. 0 means no change
	 * @return The current generation after we are finished
	 */
	uint AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f);
	uint AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f, uint threadCount);
//	uint AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f, const char *filename);
	/**
	 * @brief Empties the pool and resets the generation count
	 */
    void Clear();

	Chromosome GetExpression(uint chrIDX);

	/**
	 * @Brief Initial attempt at a single locus LOD calculation
	 * @return scores[n] = lod where n is each locus (n==diseaseLocus is 0.0)
	 */
	void CalculateLOD(vector<double>& scores, int diseaseLocus);

	/**
	 * @brief Allows an individual to be redrawn
	 * @note Only use this when an individual is known to not exist in any of the samples. 
	 */
	void ReserveIndividual(Individual* ind);

	void ResolveGenotypes(Individual *ind);

	static float defFre1;				///<Generic value used in the initial generation of loci
	static float defFre2;				///<Generic value used in the initial generation of loci
	static float errRate;				//
	static bool randomizeAlleleFreq;	///<Set by configuration to determine how to set up the loci
	
	//Growth rate details
	static float minGrowthRate;			///<Min percentage of current population (can be +/-)
	static float maxGrowthRate;			///<Max percentage of current population (can be +/-)
	static bool allowPopulationGrowth;	///<Turn growth on or off
	static uint highBufferSize;			///<The zoomed out buffer around the target block
	static uint detailsBufferSize;		///<The smaller buffer 
	static string cssFilename;			///<The file that holds the CSS information
	static uint reportWidth;			///<The max width of the report (for images to be inline/resized)
	static uint targetPop;				///<The size of the pool where advancement will not continue
	static bool UseAdamEve; 			///<2 founder chromosomes expressed expression/2 times each
	static bool UseEden;				///<2 founders who reproduce expression times to produce initial pop
	static int FounderCount;
	static double FounderDistortion;	///<Alter each founder chromosomes by this amount
	static int MaxRepeatCount;
	static double ParentDistortion;		///<Alter each parent prior to crossing over for initial population
	static double ChildDistortion;		///<Alter each child prior to adding them to the population
	static bool WritePhasedPools;		///<Turn on/off the writing of phased pools (in addition to binary)
	static size_t PhasedPoolWriteSize;	///<How many expressions to write to the phased pool

	static void SetOffspringLimits(uint min, uint max);
	/**
	 * @brief Writes LD Data around the given block
	 */
//	uint WriteLdData(const char *filename, float thresh, HaplotypeBlock *block);
//	uint WriteLdData(const char *filename, float thresh, BlockListNode *block);


	/**
	 * @brief Control the number of blocks written to the final report
	 */
	static size_t numberOfBlocksToReport;

	static float maxLdConsideration;	
	static size_t LdSpread;				///<The number of snps on either side of the block

//	static DistMappingFn *mappingFn;	///<The mapping function to be used to convert recomb fraction to bases
	static bool writeDPrimePlots;		///<Turn on/off DPrime plots	
	static bool writeRSquaredPlots;		///<Turn on/off RSquared plots 
	static bool writeLdTextReport;		///<Turn on/off Text LD report
	/**
	 * @Brief Used to determine if we open and close the pool each time we perform an advancement. 
	 * @note this could allow us to perform whole genome stuff on a single machine
	 */
	static bool closePoolBetweenAdvances;

	static void GenerateReport(ostream& os, uint headerWidth);
	
	static bool continueRunning;
	/**
	 * @brief For testing purposes, switch between traditional and the new method for pool advancements
	 */
	static bool UseOriginalCrossing;

	/**
	 * @brief allow system to turn on/off fastLD calculations
	 * @note Fast LD is used only for scanning a growth curve for the proper state
	 */
	static bool fastLD;
	static uint maxLDIndividuals;
	static uint plotScanSize;
	static uint minOffspring;
	static uint maxOffspring;



	typedef std::vector<Chromosome> PoolType;
	typedef std::vector<MinChrom> MinPoolType;
	typedef std::map<std::string, uint> IndLookupType;

	void PopulatePool(Individual* ind);
	/**
	 * @brief Calculates the allele frequencies
	 */
	bool GetAlleleFrequencies(uint locus, float &af1, float &af2);

	/**
	 * @brief Convenience overload for complete chromosome information
 	 * @param os The stream to be written to
	 * @param f The function used to determine the distance
	 * @param minAlFreq The minimum allele frequency that will be written
	 * @note This just iterates from locus 0 to locus lociCount
	 */
	uint WriteMarkerInfo(std::ostream& os, float minAlFreq=0.0);
	void WriteMendelFormat(ostream& os);
	/**
	 * @brief Write the marker information file to a stream
	 * @param os The stream to be written to
	 * @param first The first locus to be written
	 * @param last The last locus to be written
	 * @param f The Mapping function to be used. If NULL, we will use the local fn
	 * @param minAlFreq The min frequency that will be written
	 * @note Basically, work through the loci from first to last. If it is important
	 * f is NULL, the local function is used and the position will be reset to 0. 
	 */
	uint WriteMarkerInfo(std::ostream& os, uint first, uint last, float minAlFreq=0.0);

	/**
	 * @brief Write the pool in binary format to the specified file
	 * @Note Below is the file format
	 *	Binary Format: 
	 *	UINT (4 bytes) Chromosome ID
	 *	UINT (4 bytes) generation 
	 * 	UINT (4 bytes) expression count (pool size)
	 *	UINT (4 bytes) Number of loci
	 *	BLOB (determined by 2*num loci rounded up to nearest 32 bit value)
	 */
	void WriteBinaryPhased(const char *filename);

	/**
	 * @brief Write the pool in binary format to stream
	 */
	void WriteBinaryPhased(ofstream &file);

	/**
	 * @brief loads the data from binary file
	 */
	void ReadBinaryPhased(const char *filename, bool forceLoad);	//= false);
	void ReadBinaryPhased(ifstream &file);

	void ClosePool(bool keepMinimalGenotypes);

	bool LoadArchive(const char *filename);
	bool SortLoci();
	static bool UseAltLD;

	/**
	 * @brief Pass genotypes to children. Ensure that the chromosome has been fully loaded
	 */
//	void ResolveGenotypes(Sample *dataset);
	uint GetGenerationCount() { return generationCount; }
	void AddModelLocus(Locus &l);	
	void SetModelLoci(vector<uint> *ml);
	void ReadMinimalBinary(const char *filename);
	void ResetModelLoci();

	ChromosomeType GetType() { return type; }

	ChromPool::BlockDefinition &GetDefaultBlock();
	ChromPool::BlockDefinition &GetBlock(size_t idx);
	size_t GetBlockCount();

	void WriteConfiguration(ofstream &ofile);
	void ClearBlocks();
	void ClearLoci();
	double XOLambda();
private:
	bool reloadLocusFile;
	void LoadSeedData(uint expressionCount);
	bool ContinueGrowing();
	static void *PopulatePool(void *args);
	void PopulatePool(PoolType &source, PoolType &dest);

	void InitSummaryReport(const char *filename);
	bool frequenciesDetermined;
	///Used to determine the number of columns are in a line
//	uint CountColumns(const char *line);				
	

	///<Simple cheat to check if we've created this individual before
	IndLookupType usedIndividuals;	
	ChromosomeType type;				///<Distinguish between filebased and block based chromosomes (for gui)
    PoolType pool;						///<The actual pool
	MinPoolType minPool;				///<The pool of minimum alleles
	MinChrom::Common poolDetails;	///<to avoid repeating data
	uint chromID;						///<Indicate which chromosome the block is associated with
  	uint generationCount;				///<How many generations lead us to this point
	uint lociCount;						///<Number of loci associated with the local pool

	uint blockCount;					///<Used to initialize the chromosome's loci

	LocusArray loci;					///<Description of the various loci associated with the chromosomes
	vector<BlockDefinition> blockPrototypes;
	BlockDefinition &DrawBlockDefinition();

	string locSource;					///<Used for a chromosome that is to be based on a source file
	string gridSource;					///<Association grid
	string seedSource;					///<Used to indicate where seed data lives (if chromosome is seeded)
	bool positionsDefined;				///<Record that we have real positions defined
	Utility::Random rnd;				///<Random number generator
	
//	RecDistType recombinationDistribution;
	RecDistType recombIndexLookup;		///<Each index represents the locus prior to the cross over

	double poissonLambda;				///<Mean # of xo events over length of the chromosome

	ofstream *summaryReport;

	vector<uint> modelLoci;				///<Loci associated with a model on this chromosome
	bool isMinimal;
	string label;						///<This way the chromosome will be recognized
	static pthread_mutex_t poolLock;
	static pthread_mutex_t reportLock;
	int servedPools;
	string seedFmt;
	uint headerCount;

};




class PoolArchiveException : public Utility::Exception::Base {

	PoolArchiveException(const char *filename) : filename(filename) { } 
	PoolArchiveException(const string& filename) : filename(filename) { }
	PoolArchiveException(const PoolArchiveException& other) : filename(other.filename) { }
	virtual ~PoolArchiveException() { cout<<"PoolArchiveException("<<filename<<")\n";}

	virtual string GetErrorMessage() const;
protected:
	string filename;
};

inline
bool ChromPool::ContinueGrowing() { 
	bool doContinue = false;
	LOCKPOOL;
	if (servedPools-- > 0)
		doContinue = true;
	UNLOCKPOOL;
	return doContinue;
}

inline
string PoolArchiveException::GetErrorMessage() const {
	return "A problem was encountered with the archive, "+string(filename);
}

inline
uint ChromPool::BlockDefinition::GetSnpCount() {
	//Figure out how many loci there are locally
	int diff = maxSnpCount - minSnpCount;
	return minSnpCount + Utility::Random::globalGenerator(diff);
}

inline
float ChromPool::BlockDefinition::GetBlockMapDistance() {
	//First one we use the block strength variation from 0.5
	float rDiff = maxBlckMap - minBlckMap;
	return minBlckMap + (float)(Utility::Random::globalGenerator.drand() * rDiff) + 0.0000001;
}

inline
float ChromPool::BlockDefinition::GetSnpMapDistance() {
	float rDiff = maxSnpMap - minSnpMap;
	return minSnpMap + (float)(Utility::Random::globalGenerator.drand() * rDiff)+ 0.0000001;
}

inline
void ChromPool::BlockDefinition::WriteConfiguration(ostream &str) {
	str<<"ADD_BLOCK "<<minSnpCount<<" "<<maxSnpCount<<" "<<minBlckMap<<" "<<maxBlckMap<<" "<<minSnpMap<<" "<<maxSnpMap<<" "<<frequency<<"\n";
}


}

#endif
