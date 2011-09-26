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
#include "diseasemodel.h"
#include "locreport.h"

#ifdef USE_XY
#include "cpoolxy.h"
#endif //USE_XY
namespace Simulation {

class Sample;
class StatusModel::DiseaseModel;
using namespace Simulation::Visualization;

#define POOLLOCK 	pthread_mutex_lock(&chromPoolLock)
#define POOLUNLOCK 	pthread_mutex_unlock(&chromPoolLock)

#define SUMMARYLOCK pthread_mutex_lock(&summaryLock)
#define SUMMARYUNLOCK pthread_mutex_unlock(&summaryLock)


typedef vector<Simulation::Visualization::LocusSelection > LocusSelectionArray;
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
	bool LoadPoolContents(const char *project, uint generation, bool closeAfterLoading);

	bool LoadPoolContents(const char *project, uint generation, uint chrID, bool closeAfterLoading);

	bool LoadLoci(const char *project, uint generation);
	bool LoadLoci(const char *project, uint generation, uint chrID);

	bool Close();
	/**
	 * @brief Write the pools to their respective file
	 */
	bool Dump(const char *prj);


	string ProduceOverviewLD(const char *prj, bool forDataset, bool loadFirst, ofstream &summary, const char *growthChart);

	string GenerateSampleMarkerInfo(const char *prj, float minFreq);

	/**
	 * @brief Generate the filename associated with the whole chromosome marker information
	 * @param pos The chromosome ID
	 * @param prj The project name used as part of the filename
	 */
	string GenerateMarkerInfoFilename(uint pos, uint chr, const char *prj, bool forDataset);

	/**
	 * @brief Generate the filename associated with a subset from an overside pool
	 * @param pos The Chromosome ID
	 * @param segmentID Which piece this file represents
	 * @param prj The project name (used as the first portion of the filename)
	 */
	string GenerateMarkerInfoFilename(uint pos, uint chr, const char *prj, uint startLoci, bool forDataset);

	/**
	 * @brief Add a pool to the manager
	 * @note This puts a pool completely under local management- including delete
	 */
	void AddPool(ChromPool *pool);

	ChromPool *GetChromosome(uint idx);
	void DeleteChromosome(uint idx);
	
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
	bool AdvanceGenerations(uint moreGenerations, PopulationGrowth::GrowthRate *f, bool loadFirst);
	
	/**
 	 * @brief Iterate through pools and calculates each pool's allele freqs
	 */
	void CalculateAlleleFrequencies();
	
	/**
	 * @brief Create a new chromosome pool to be managed locally
	 */
	ChromPool *AddChromosome(uint blockCount, ChromPool::BlockDefinition &defaultBlock, const char *label);
	ChromPool *SeedChromosome(const char *seed, const char *loc, const char *label, const char *fmt, uint headerCount);
	ChromPool *AddChromosome(const char *locFilename, const char *label);
	ChromPool *AddChromosome(const char *locfilena, const char *gridFilename, const char *label);

#ifdef USE_XY
	void AddChromosomeXY(const char *locFilename, const char *label, Utility::Random& rnd);	
	AlleleSource<LocusXY> *DrawX();
	AlleleSource<LocusXY> *DrawY();
#endif //USE_XY
	/**
	 * @brief Brings the pools up with the proper data according to the generation
	 */
	bool InitializePools(uint startGen, uint poolSize, const char *prj, bool doLoad, bool doDump, bool closeAfterCreating);
	bool InitializeLoci(uint startGen, const char *prj, bool doLoad);

	/**
	 * @brief Forces the allele frequencies associated with a given chromosome and locus
	 */
	bool ForceAlleleFrequency(uint chrID, uint locusID, float al1, float al2);

	bool ForceAlleleFreqRange(uint chrID, uint locusID, float al1, float al2);

	/**
	 * @brief Seek out the frequency for a given locus on a given chromosome and return it (by ref)
	 */
	bool GetAlleleFrequency(uint chrID, uint locId, float &al1, float &al2);


	LocusMap &GetLocusMap() { return locusMap; }

	/**
	 * @brief Start with a random population. This is required if you don't load one from scratch
	 * @param populationCount The number of instances to be randomly created
	 * @param doDump Write the population immediately upon creation
	 * @param closeAfterCreating Archive to file and unload afterward? This overrides the doDump
	 */
	void CreateInitialPopulation(uint populationCount, bool doDump, bool closeAfterCreating);


	void DrawIndividual(Individual& ind, bool isXX);
	/**
	 * @brief Randomly apply gender (if XY exists)
	 */
	void DrawIndividual(Individual& ind);
	/**
	 * @brief Returns the number of chromosome pools are under management
	 */
	size_t Size() { return pools.size(); }

	//void WriteLocusFile(std::ostream& os, uint first, uint last, ChromPool *ch);
	
	/**
 	 * @brief Generate the phased output from the local set of chromosomes. 
 	 * @param filename Base filename to be used as first portion of the files
	 * @param produceWindows If true, detailed windows of the chromosome's pieces will be
	 						 written.
	 */
//	bool WritePhasedData(const char *filename, bool produceWindows, bool forDataset);


	/**
	 * @brief Generate phased output from each of the chromosomes as one file/chromosome
	 */
	string WritePhasedOverview(const char *project, bool forDataset);
	bool WritePhasedWindowed(const char *project, bool forDataset);
//	bool WriteLdInputForHaploview(const char *project, bool loadFirst);
	string WriteMarkerInfoOverview(const char *project, bool forDataset);
//	bool WriteMarkerInfoWindowed(const char *project, bool forDataset);

	//bool ProduceOverviewLD(const char *phData, const char *markerInfo, bool forDataset);
//	bool ProduceWindowedLD(const char *project, bool forDataset);
	

	bool GetLocusReference(const char *label, Locus& locus);

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
	
	static bool generateOverviewLD;

	//static bool useDistanceMaps;
	//static void* distanceMappingFunction;
	

	size_t GetPoolCount();

	size_t GetExpressionCount();

	void InsertModel(StatusModel::DiseaseModel *model);

	/**
	 * @brief Perform necessary loading/unloading required for dataset generation
	 * @note Basically, this will load required chromosomes into memory 
	 */
	bool PrepForSampling(const char *project, uint generation, StatusModel::DiseaseModel *model);

	bool ResolveGenotypes(const char *project, uint generation, Sample *samplePopulation, StatusModel::DiseaseModel *model);
	bool ResolveGenotypes(vector<Individual*>& people, StatusModel::DiseaseModel *model);
	void SetLocusSelectors(LocusSelectionArray& locSel);

	static uint maxThreadCount;
	ChromPool *GetNextPool(uint &idx);
	static void *ProduceOverviewLD(void *arg);
	static void *AdvanceGenerations(void* args);

	static uint simultaneousChrom;
	static uint threadsPerChrom;
	void SetProjectName(const char *projectName);
protected:

	bool VerifyHaploview();

	void SetPoolGeneration(uint generation);	


	/**
	 * @brief Produce a filename based on a few details of the state of the application or the pools themselves
	 */
	string GeneratePhasedFilename(uint gen, uint position, const char *prj, bool forDataset);

	string GeneratePhasedFilename(uint gen, uint pos, const char *prj, int start, bool forDataset);

	vector<ChromPool *> pools;					///<The various chromosome pool objects being managed

#ifdef USE_XY
	CPoolXY *poolxy;
	LocusManager<LocusXY> *locusXY;				///<XY Loci
#endif //USE_XY
	uint currGeneration;						///<Which generation we are currently at
	bool poolIsPopulated;						///<Quick test to avoid running on empty pools
	LocusMap locusMap;							///<Used to map SNP Labels to chromosome/locus 

	LocusReport locusReport;	
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
	string projectName;							///<Used as part of filenames


	
	uint threadedPoolIdx;						///<This will be used for navigating pools during threaded efforts
	static pthread_mutex_t chromPoolLock;		///<Lock down the pool whenever we are doing threaded stuff with it
	static pthread_mutex_t summaryLock;			///<Lock down the summary stream

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
