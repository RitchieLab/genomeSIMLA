//
// C++ Interface: poolmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPOOLMANAGER_H
#define SIMULATIONPOOLMANAGER_H

#include "utility/types.h"
#include "chrompool.h"
#include "growthrate.h"

namespace Simulation {

class Sample;

/**
@brief A helper class to manage the different chromosome pools
Use this class to setup each of the chromosomes or the respective proxies assocaited with remote chromsome pools

Most of the control required of a specific pool can be applied to the manager object- and it will make the necessary
calls to each of the children pools

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PoolManager{
public:



    PoolManager();

    ~PoolManager();

	/**
	 * @brief Have each pool load its data from a file. 
	 * The filename is correctly generated according to the project, generation # and pool details. 
	 * It is assumed that pools will be loaded in the same order as the were created on previous runs
	 * in case that position is used as part of the filename
	 * @return true/false indicating successful load. This should represent a high confidence that
	 * all data was succesfully loaded, even if the pools reside on remote machines. If the remote
	 * machine never responds, it is possible we will time out locally...but that has yet to be seen
	 */
	bool Load(const char *project, uint generation);

	/**
	 * @brief Write the pools to their respective file
	 */
	bool Dump(const char *prj);

	bool GenerateLDMaps(const char *prj, bool leaveLocusFiles);
	bool GenerateLDMap(const char *settings, const char *filename, const char *markerInfo, bool leavePhased);
	bool GenerateLDMap(const char *filename, const char *markerInfo, bool leavePhased);
	bool GenerateLDMapOverview(const char *filename, const char *markerInfo, bool leavePhased);

	/**
	 * @brief Generate the filename associated with the whole chromosome marker information
	 * @param pos The chromosome ID
	 * @param prj The project name used as part of the filename
	 */
	string GenerateMarkerInfoFilename(uint pos, const char *prj);

	/**
	 * @brief Generate the filename associated with a subset from an overside pool
	 * @param pos The Chromosome ID
	 * @param segmentID Which piece this file represents
	 * @param prj The project name (used as the first portion of the filename)
	 */
	string GenerateMarkerInfoFilename(uint pos, uint segmentID, const char *prj);

	/**
	 * @brief Add a pool to the manager
	 * @note This puts a pool completely under local management- including delete
	 */
	void AddPool(ChromPool *pool);

	/**
	 * @brief Returns the current generation
	 */
	uint GetCurrentGeneration();

	/**
	 * @brief Sends configuration information about the chromosomes associated with the manager to the stream
	 */
	void GenerateReport(ostream &os, uint headerWidth);
	
	/**
	 * @brief simulate the progress of moreGenerations generations
	 */
	bool AdvanceGenerations(uint moreGenerations, PopulationGrowth::GrowthRate *f);

	/**
	 * @brief Create a new chromosome pool to be managed locally
	 */
	ChromPool *AddChromosome(uint blockCount, float minR, float maxR, ChromPool::BlockDefinition &defaultBlock);

	ChromPool *AddChromosome(const char *locFilename);

	/**
	 * @brief Brings the pools up with the proper data according to the generation
	 */
	bool InitializePools(uint startGen, uint poolSize, const char *prj);

	/**
	 * @brief Forces the allele frequencies associated with a given chromosome and locus
	 */
	bool ForceAlleleFrequency(uint chrID, uint locusID, float al1, float al2);

	bool ForceAlleleFreqRange(uint chrID, uint locusID, float al1, float al2);

	/**
	 * @brief Seek out the frequency for a given locus on a given chromosome and return it (by ref)
	 */
	bool GetAlleleFrequency(uint chrID, uint locId, float &al1, float &al2);
	/**
	 * @brief Start with a random population. This is required if you don't load one from scratch
	 */
	void CreateInitialPopulation(uint populationCount);

	/**
	 * @brief Returns the number of chromosome pools are under management
	 */
	uint Size() { return pools.size(); }

	//void WriteLocusFile(std::ostream& os, uint first, uint last, ChromPool *ch);
	
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
		ChromPool *GetNext();			
		~Iterator();									///<Destruction
		void Reset();									///<Reset the iterator
		Iterator(vector<ChromPool *> *repository);		///<Basic contruction
	protected:	
		vector<ChromPool *> *repository;				///<The repository we are iterating through
		vector<ChromPool *>::iterator position;			///<The current position
		
	};


	/**
	 * @brief Acquire the iterator for the first Chromosome Pool
	 * @note This iterator doesn't guarantee that it will traverse the chromosomes in the same order each time 
	 */
	Iterator GetIterator() {
		return Iterator(&pools);
	}
	
	static string pathToJava;
	static string pathToHaploview;
	static string haploviewSettings;	
	static string haploviewOverviewSettings;
	static string javaSettings;							//This can be used to change the memory assignment

	static int haploWindowSize;
	static int haploWindowStride;
	static bool generateOverviewLD;

	//static bool useDistanceMaps;
	//static void* distanceMappingFunction;
	

protected:

	bool VerifyHaploview();


	/**
	 * @brief Produce a filename based on a few details of the state of the application or the pools themselves
	 */
	string GeneratePhasedFilename(uint gen, uint position, const char *prj);

	string GeneratePhasedFilename(uint gen, uint pos, const char *prj, int start);

	/**
	 * @brief Produce a filename that is appropriate for the locus report. 
	 * This filename will be unique for generation and position (based on the project name)
	 */
	string GenerateLocusFilename(uint gen, uint position, const char *prj);

	vector<ChromPool *> pools;					///<The various chromosome pool objects being managed
	uint currGeneration;						///<Which generation we are currently at
	bool poolIsPopulated;						///<Quick test to avoid running on empty pools

	/**
	 * @brief Store the details required for explicit allele frequency definition
	 */
	struct ForcedAF {
		uint chrID;								///<ID of the chromosome (0 based)
		uint locus;								///<Index of the locus (0 based)
		float af1;								///<frequency for 1rst allele
		float af2;								///<frequency for the 2nd allele
		ForcedAF() : chrID(0), locus(0), af1(0.0), af2(0.0) {}
		ForcedAF(uint c, uint l, float af1, float af2) : chrID(c), locus(l), af1(af1), af2(af2) { }
	};

	vector<ForcedAF> forcedFrequencies;			///<Cache the forced frequencies

};

inline
uint PoolManager::GetCurrentGeneration() { return currGeneration; }




inline
void PoolManager::Iterator::Reset() {
	position = repository->begin();
}

inline
PoolManager::Iterator::Iterator(vector<ChromPool *> *repository) : repository(repository) {
	position = repository->begin();
}



inline
PoolManager::Iterator::~Iterator() { }

inline
ChromPool *PoolManager::Iterator::GetNext() {
	ChromPool *current = NULL;
	if (position != repository->end()) { 
		current=*position;
		position++;
	}
	return current;
}


}	



#endif
