//
// C++ Interface: cpool
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONCPOOL_H
#define SIMULATIONCPOOL_H
#include "utility/types.h"
#include "growthrate.h"
#include "locusmanager.h"
#include "utility/random.h"
#include "allelesource.h"
#include "utility/stat.h"
#include <pthread.h>
#include "locreport.h"
#include "utility/array2d.h"

#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif

namespace Simulation {
	
template <class T>
struct Population {
	typedef vector<AlleleSource<T> *> Type;
};


#define LOCKPOOL 	pthread_mutex_lock(&poolLock)
#define UNLOCKPOOL 	pthread_mutex_unlock(&poolLock)

#define LOCKREPORT  pthread_mutex_lock(&reportLock)
#define UNLOCKREPORT pthread_mutex_unlock(&reportLock)
	
extern unsigned int targetPop;					///<Target population 
extern unsigned int numberOfBlocksToReport;		///<Control the number of blocks written to the final report
extern unsigned int highBufferSize;				///<The zoomed out buffer around the target block
extern unsigned int detailsBufferSize;			///<The smaller buffer 


/**
@brief New Chromosome pool, using Allele sources instead of chromosomes. 
This new version should work more efficiently with generational advancements and could potentially represent a significant savings of disk-space, depending on the length of simulations

Changes from the original chromosome:
	Block based chromosomes don't mean anything!		-- Well, there should be a specialized BB conversion system that
		creates locus arrays for chromosomes. The chromosome pool shouldn't be a swiss army class. 
	XO events are recorded for every generation and chromosomes must be "realized" prior to writing datasets. This means that
		increasing populations won't be a linear growth, but there will be logarithmic growth that eventually will surpass the
		linear growth that the previous mechanism used.The "Source" population must be exist in RAM when writing genotypes
		to file
	Pools must be responsible for "Minimizing" a given allele source, if we don't want the entire contents to remain in memory. 
		The previous structure knew how to record only model loci, but this was changed for simplifying the code. 
		Also, the source generation will be part of the binary format. This pool must have valid pointers at the time of 
		advancement, though, we will be able to load individual pools after cross over at any point prior to realization in 
		order to realize chromosomes. In the case of generational advancement, this will never be necessary unless we need 
		to "Bake" the pools in order to reset the memory back to the linear growth point where the new method exceeds 
		the previous method and the occassional analysis (LD). To resume a minimal memory footprint and to resume a 
		speedy advancement rate, we could "bake" the contents at an analysis dump, and reset those XO structures back to
		nothing...this would cause a loss in allele tracking data, but would keep the advancement rate to reasonably high
		rates...possibly high enough rates to realistically handle mutation without cluster resources.
	Loci are managed by a Locus Manager object.
		This goes back to the BB chromosomes above, and also handling tasks like loading locus information from file. This will be important in handling the XY oddities without disrupting the chromosome too much
	Analysis is now somewhat trickier
		Since raw data exists only when the source pool is fully in memory, things like allele frequencies can't be evaluated until we stop and reload the source pool(s). This is time consuming and should only be performed when necessary. 
	Certain aspects of the pool will stay in memory until certain functions are called.
		These include project details (for filenames) as well as locus details. This will prevent a lot of the reloading and weird behavior that haunts the old chromosome pools. So, unless a Reset function is called, loci will not be "Read" but "reloaded". 
	Simplified functions
		Instead of providing an "adapter like" approach for accessor functions related to the pool data, the data will be available via a single function call. This is somewhat dangerous, but should reduce the volumes of dumb accessor functions. An option would be to use a pool container with two levels of accessor functions. We could expose a "light" version to protect the data somehow....need to think about that.

	@NOTE During "Suspension", the source population is the only thing that is cleared....so, as we approach later generations, the memory footprint will continue to grow. If this is a problem, we might have to write these inherited alleleSources to file and reload them so that we can work within memory constraints
	@author Eric Torstenson
*/
template <class T>
class CPool {
public:
	/**
	 * @brief Constructor...
	 * @param chromID Index of this chromosome
	 * @param project Prefix used to write/read all files
	 * @param label Used to help user identify a given chromosome in reports and filenames
	 */
    CPool(int chromID, const char *project, const char *label);

	/**
	 * 
	 */
    ~CPool();

	/**
	 * @Brief Advance by population by specified number of generations
	 * @param generations Number of generations to advance
	 * @param f Growth function (determine what size to grow to
	 * @param tCount number of threads to run in
	 */
	uint AdvanceGenerations(uint generations, PopulationGrowth::GrowthRate* f, Utility::Random& rnd, uint tCount);

	/**
	 * @brief Create new chromosomes, based on the locus array
	 */
	virtual void Init(LocusManager<T>* loci, Utility::Random& rnd);
	
	void Open(LocusManager<T>* loci, Utility::Random& rnd);			///<We should always start by loading generation 0 (source population) if such a pool exists

	void RefreshSources();
	void Refresh(int generation);			///<Load the current generation
	
//	void Close();							///<Close pools entirely (freeing all memory) including the source population
	
	void Save(int generation = -1);			///<Write current population to pool

	/**
	 * @Brief Override the project name with a new setting
	 */
	void SetProject(const char *project); 

	void Bake();							///<Realizes current population

	void Wake();

	void Suspend();							///<Mimize pool memory by removing genotype data
	void Suspend(typename Population<T>::Type *pop);	
	/**
	 * @Brief Create initial population using random draws based on minor allele frequency
	 */
	virtual void BuildInitialPopulation(Utility::Random& rnd, int expressionCount);

	/**
	 * @brief Populate pool based on source population
	 */
	virtual void PopulatePool(typename Population<T>::Type* source, typename Population<T>::Type* dest, Utility::Random& rnd);

	/**
	 * @Brief Threaded pool population (generational advancement)
	 */
	static void *PopulatePool(void *args);

	bool ContinueGrowing( );

	AlleleSource<T> *At(int idx);

	pthread_mutex_t poolLock;				///<Multithread lock for pool operations	
	pthread_mutex_t reportLock;		///<Multithread lock for reports		

	virtual int GetPopulationSize();

	std::string GetLabel();					///<Returns descriptive label associated with a given pool

	virtual void WriteLdData(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh );

	size_t CalculateAlleleFrequencies(uint maxIndividualCount);

	int GetCurrentGeneration(); 

	bool GetAlleleFrequencies(uint locus, float& af1, float& af2);

	void SetModelLoci(vector<uint>& modelLoci);

	AlleleSource<T> *Draw(Utility::Random& rnd);				///<Draw a chromosome from the pool (for now, not crossed)

	virtual void DebugPrint(ostream& os, int start, int stop);

protected:
	/**
	 * @brief Release memory associated with a population
	 */
	void Purge(typename Population<T>::Type *population);
	string GetFilename(int generation);

	/**
	 * @Brief fill events with cross over indices usin rnd as random number generator, based on poisson lambda
	 */
	bool BuildXOEvents(Utility::vector<size_t>& events, Utility::Random& rnd, float poisson=0.0);
	bool BuildXOEvents(Utility::vector<size_t>& events, Utility::Random& rnd, bool isMale);	///<Builds XO points where lambda is adjusted for gender 

	int curGeneration;						///<Current generation
	bool isSuspended;						///<Indicate that pointers are valid, but no genotypes are avaialble (except model related)
	int id;									///<Chromosome index
	std::string project;					///<cached project
	std::string label;						///<Id written to locus files
	/**
	 * @brief These represent the "founders"..and potentially contain complete genotype information
	 */
	typename Population<T>::Type *sourcePopulation;
	
	Utility::RBTree<size_t, AlleleSource<T>*> sourceAlleles;	
	
	/**
	 * @brief Population at current generation. If This is empty, the source population IS the current generation
	 */
	typename Population<T>::Type *curPopulation;
	LocusManager<T> *loci;					///<The locus array associated with this chromosome
	float poissonLambda;					///<Baseline lambda for chromosome. Gendered XO will adjust this accordingly
	float gdStart;							///<The genetic distance between the beginning chromosome and the first SNP
	float gdLength;							///<Length of the chromosome in map units
	int remainingChroms;					///<Number of chromosomes to be added to a growing pool
 	uint generationCount;					///<How many generations lead us to this point
	vector<uint> modelLoci;					///<The disease loci associated with this chromosome	
	int sourceGeneration;					///<This indicates which generation holds that raw data necessary to realize the current population
};

template <class T>
inline
void CPool<T>::SetModelLoci(vector<uint>& modelLoci) {
	this->modelLoci = modelLoci;
}

template <class T>
inline
bool CPool<T>::GetAlleleFrequencies(uint locusIdx, float& af1, float& af2) {
	typename Population<T>::Type *pop = curPopulation;
	if (curPopulation == NULL) 
		pop = sourcePopulation;

	//For now, this is only used for disease models. As a result, we want to take the average
	//of the two frequencies.
	T locus;
	loci->At((int)locusIdx)->Distill(locus);
	af1 = locus.Freq1();
	af2 = locus.Freq2();	
	return true;
}

template <class T>
inline
void CPool<T>::SetProject( const char *project) {
	this->project =project;
}

template <class T>
inline
int CPool<T>::GetPopulationSize() {
	if (curPopulation)
		return curPopulation->size();
	else if (sourcePopulation) 
		return sourcePopulation->size();
	return 0;
}

template <class T>
inline
AlleleSource<T> *CPool<T>::At(int idx) {
	if (curPopulation) {
		assert(idx < curPopulation->size());
		return curPopulation->at(idx);
	}
	else {
AlleleSource<T> *crom = sourcePopulation->at(idx);
		assert(idx < sourcePopulation->size());
		return sourcePopulation->at(idx);
	}
}	

/*****/
template<class T>
inline
string CPool<T>::GetFilename(int generation) {
	stringstream ss;
	ss<<project<<"."<<label<<"."<<generation<<".phased";
	return ss.str();
}

template <class T>
inline
void CPool<T>::DebugPrint(ostream& os, int start, int stop) {
	typename Population<T>::Type *pop = curPopulation;
	if (curPopulation == NULL) 
		pop = sourcePopulation;
	
	for (int i=start; i<stop; i++) {	
		os<<"-- ";
		pop->at(i)->ShowGenotypes(os, 10, '\t');
		os<<"\n";
	}
}

template <class T>
inline
int CPool<T>::GetCurrentGeneration() { 
	return curGeneration;
}


template<class T>
inline
void CPool<T>::Open(LocusManager<T> *loci, Utility::Random& rnd) {
	Init(loci, rnd);
	if (sourcePopulation) 
		Purge(sourcePopulation);
	else
		sourcePopulation = new vector<AlleleSource<T> *>();
	
	ifstream file(GetFilename(0).c_str());
	file.read((char*)&id, 4);
	int gen = 0;
	file.read((char*)&gen, 4);

	curGeneration = gen;
	assert(curGeneration == 0);
	size_t popSize = 0;
	file.read((char*)&popSize, 4);
	int locCount = 0;
	file.read((char*)&locCount, 4);
	assert(locCount == loci->LocusCount());
	int popType = 0;
	file.read((char*)&popType, 4);
	assert(popType == 0);
	for (int i=0; i<popSize; i++) {
		AlleleSourceSingle<T> *chrom = new AlleleSourceSingle<T>(loci, i);
		chrom->ReadBinary(file, NULL);		///<Source alleles is irrelavent for this pool 
		sourcePopulation->push_back(chrom);
	} 
}


template <class T>
inline
void CPool<T>::RefreshSources() {
	ifstream file(GetFilename(0).c_str());
	file.read((char*)&id, 4);
	int curGeneration = 0;
	file.read((char*)&curGeneration, 4);
	assert(curGeneration == 0);
	size_t popSize = 0;
	file.read((char*)&popSize, 4);
	int locCount = 0;
	file.read((char*)&locCount, 4);
	assert(locCount == loci->LocusCount());
	int popType = 0;
	file.read((char*)&popType, 4);

	for (int i=0; i<popSize; i++) {
		AlleleSource<T> *chrom = sourcePopulation->at(i);
		///<Source alleles is irrelavent for this pool 
		chrom->ReadBinary(file, &sourceAlleles);		
	} 
}

template <class T>
inline
AlleleSource<T> *CPool<T>::Draw(Utility::Random& rnd) {
	typename Population<T>::Type *pop = curPopulation;
	if (pop == NULL) 
		pop = sourcePopulation;

	AlleleSource<T> *source = pop->at(rnd((long int)pop->size()));
	//we want to create a basic copy, with no allelic stuff to copy, so that nothing is wasted if we can't use it
	vector<size_t> empty;
	return source->Cross(source, empty);
}

template <class T>
inline
std::string CPool<T>::GetLabel() {
	return label;
}


template <class T>
inline
void CPool<T>::Refresh(int generation) {
	if (curPopulation) 
		Purge(curPopulation);
	else 
		curPopulation = new vector<AlleleSource<T> *>();
	ifstream file(GetFilename(generation).c_str(), ios::binary);
	if (file.good()) {
		file.read((char*)&id, 4);
		file.read((char*)&curGeneration, 4);
		assert(curGeneration == generation);
		size_t popSize = 0;
		file.read((char*)&popSize, 4);
		int locCount = 0;
		file.read((char*)&locCount, 4);
		assert(locCount == loci->LocusCount());
		int popType = 0;
		file.read((char*)&popType, 4);
	
		for (int i=0; i<popSize; i++) {
			AlleleSource<T> *chrom;
			if (popType == 0)
				chrom = new AlleleSourceSingle<T>(loci, i);
			else 
				chrom = new AlleleSourceInh<T>(loci, i);
			///<Source alleles is irrelavent for this pool 
			chrom->ReadBinary(file, &sourceAlleles);		
			curPopulation->push_back(chrom);
		} 
	
		if (DoSuspend)
			Suspend(curPopulation);
	}
	else {
		stringstream ss;
		ss<<"Error when trying to refresh chromosome pool, "<<label<<". The file "<<GetFilename(generation)<<" doesn't seem to be valid.";
		throw Utility::Exception::General(ss.str().c_str());
	}
}	

template<class T>
inline
void CPool<T>::Save(int generation) {
	if (generation == -1)
		generation = curGeneration;
	ofstream file(GetFilename(generation).c_str(), ios::trunc|ios::binary);
	//Write all the header details
	file.write((char*)&id, 4);
	file.write((char*)&generation, 4);
	size_t popSize = 0;
	typename Population<T>::Type *population;
	if (generation == sourceGeneration) {
		population = sourcePopulation;
		assert(curPopulation==NULL);
	} else {
		population = curPopulation;
	}
	
	popSize = population->size();

	file.write((char*)&popSize, 4);
	assert(loci->LocusCount() > 0);
	int locCount = loci->LocusCount();
	file.write((char*)&locCount, 4);
	int popType = generation > sourceGeneration;
	file.write((char*)&popType, 4);
	//We only really want to save the "current" generation, whichever that is.
	//Copy this stuff from the original chrompool object
	typename Population<T>::Type::iterator itr = population->begin();
	typename Population<T>::Type::iterator end = population->end();
	while (itr != end) {
		(*itr)->WriteBinary(file);
		itr++;
	}

	LocusManagerFileBased<T> *l = (LocusManagerFileBased<T>*)loci;
	l->Save(generation, label.c_str());

}



template <class T>
inline
size_t CPool<T>::CalculateAlleleFrequencies(uint maxIndividualCount) {
	
	typename Population<T>::Type *pop = curPopulation;
	if (pop == NULL)
		pop = sourcePopulation;

	uint indCount = pop->size();
	size_t lociCount = loci->LocusCount();

	//If we want to perform "sampled" frequencies, we do that here
	typename Population<T>::Type::iterator start = pop->begin();
	typename Population<T>::Type::iterator end   = start;

	if (indCount == 0)
		return 0;
	if (maxIndividualCount == 0 || maxIndividualCount > indCount) {
		maxIndividualCount = indCount;	
		end = pop->end();
	}
	else {
		for (int i=0; i<maxIndividualCount; i++) 
			end++;
	}

	
	uint fixedAlleles = 0;
	
	Array2D<int> freq(lociCount, 2);
	typename Population<T>::Type::iterator chrom = start;
	while (chrom != end) {
		for (uint i=0; i<lociCount; i++) 
			freq(i, chrom->At(i))++;
		chrom++;
	}

	for (int i=0; i<lociCount; i++) {
		float freq1 = freq(i,0);
		if (freq1 == 1.0 || freq1 == 0.0) 
			fixedAlleles++;
		double al1 = freq1/(double)indCount;
		double al2 = freq(i, 1)/(double)indCount;
		Locus *locus = loci->At(i);
		locus->AssignFreq(al1, al2);
	}

	if (fixedAlleles == lociCount) {
		chrom = start;
		while (start != end){
			cerr<<"Locus Count ("<<lociCount<<")\t";
			chrom->ShowGenotypes(cerr, 20, '>');
			cerr<<"\n";
		}
		cerr<<"We have a problem. There are no variants in our pool at all. Every one of our "<<indCount<<" expressions is the same\n";
		abort();
	}

	return fixedAlleles;
}


template<class T>
inline
void CPool<T>::Wake() {
	if (DoSuspend) {
		RefreshSources();
		Refresh(curGeneration);
	}
cerr<<"CPool::Wake\n";
}

template<class T>
inline
void CPool<T>::Suspend(typename Population<T>::Type *pop) {
	typename Population<T>::Type::iterator itr = pop->begin();
	typename Population<T>::Type::iterator end = pop->end();
	while (itr!=end) {
		(*itr)->Suspend(modelLoci);
		itr++;
	}
cerr<<"Cpool::Suspend\n";
}

template<class T>
inline
void CPool<T>::Suspend() {
	//For now, this is the same, but I can see it changing slightly, if we use suspension during advancement and need to use temporary files
	if (DoSuspend) {
		Save(0);
		if (curGeneration > 0)
			Save(curGeneration);
		Suspend(sourcePopulation);
		Suspend(curPopulation);
		isSuspended = true;
	}
}

/**
 * We use bake so that we have complete genotyping information at the current generation, with no dependancies on a prior generation...such as in preparation of datasets
 */
template<class T>
inline
void CPool<T>::Bake() {
	//This makes no sense, if the pool hasn't been opened
//	assert(isOpen);
	assert(sourcePopulation);
	
	int sourceID = sourceAlleles.FindMax()->GetKey();
	sourceAlleles.Clear();

	if (isSuspended) 
		Wake();
	
	if (curPopulation) {
		int count = curPopulation->size();
		
		for (int i=0; i<count; i++) {
			AlleleSource<T> *originalSource = curPopulation->at(i);
			AlleleSource<T> *realizedSource = originalSource->Realize(++sourceID);
			delete originalSource;
			curPopulation->at(i, realizedSource);
			sourcePopulation->Insert(sourceID, realizedSource);
		}
		
		//Get rid of the source population
		typename Population<T>::Type::iterator itr = curPopulation->begin();
		typename Population<T>::Type::iterator end = curPopulation->end();
		while (itr != end) {
			delete *itr;
			itr++;
		}
		Purge(sourcePopulation);
		sourcePopulation = curPopulation;
		sourceGeneration = curGeneration;
		curPopulation = NULL;
	}
	
}

template<class T>
inline
void CPool<T>::Init(LocusManager<T>* loci, Utility::Random& rnd) { 
	this->loci = loci;
	gdStart = loci->At(0)->MapPosition();
	int locusCount = loci->LocusCount();
	gdLength = loci->At(locusCount -1)->MapPosition()-gdStart;
	poissonLambda = gdLength * 0.01;
}

template<class T>
struct AdvancementArg {
	typename Population<T>::Type *sourcePop;
	typename Population<T>::Type newPop;
	Utility::Random rnd;
	uint count;
	CPool<T> *pool;

	AdvancementArg(typename Population<T>::Type *sourcePop, uint count, CPool<T> *pool, Utility::Random& rnd) : 
					sourcePop(sourcePop), rnd(rnd), count(count), pool(pool) { 
		newPop.reserve(count);
	}
};

template<class T>
inline
bool CPool<T>::ContinueGrowing( ) {
	bool doContinue = false;
	LOCKPOOL;
	doContinue = remainingChroms-- > 0;
	UNLOCKPOOL;
	return doContinue;
}
template<class T>
inline
void CPool<T>::PopulatePool(typename Population<T>::Type *source, typename Population<T>::Type *dest, Utility::Random& rnd) {
	int originalSize = source->size();

	size_t eventsObserved = 0;
	while (ContinueGrowing()) {
		int matID = rnd(originalSize);
		int patID = rnd(originalSize);
		vector<size_t> events; 
		BuildXOEvents(events, rnd);
		dest->push_back(source->at(patID)->Cross(source->at(matID), events));
	}
}



template<class T>
inline
void *CPool<T>::PopulatePool(void *args) {
	AdvancementArg<T> *thArg = (AdvancementArg<T>*)args;
	thArg->pool->PopulatePool(thArg->sourcePop, &thArg->newPop, thArg->rnd);

	return args;
}

template<class T>
inline
void CPool<T>::Purge(typename Population<T>::Type *population) {
	if (population) {
		typename Population<T>::Type::iterator itr = population->begin();
		typename Population<T>::Type::iterator end = population->end();
	
		while (itr != end) {
			delete *itr;
			itr++;
		}
		population->clear();
	}
}
template<class T>
inline
uint CPool<T>::AdvanceGenerations(uint generations, PopulationGrowth::GrowthRate* f, Utility::Random& rnd, uint tCount) {
	tCount--;
	
	bool continueRunning = true;
	uint goalPop = targetPop;
	if (goalPop == 0)
		goalPop = (uint)-1;

	typename Population<T>::Type *newPopulation = new vector<AlleleSource<T> *>();

	typename Population<T>::Type *cur = curPopulation;
	if (cur == NULL)
		cur = sourcePopulation;
		//cur = new vector<AlleleSource<T>*>(*sourcePopulation);
	uint popSize = cur->size();
	typename Population<T>::Type *next = newPopulation;
	
	for (uint currGen = 0; currGen < generations && popSize < goalPop && continueRunning; currGen++) {
		if (next == sourcePopulation) 
			next = new vector<AlleleSource<T>*>();
		else 
			Purge(next);

		uint endSize = (*f)(++generationCount);
		
		if (endSize > goalPop)
			endSize = goalPop;

		remainingChroms = endSize;
		cout<<".";cout.flush();

		pthread_t threads[tCount];
		for (uint th=0; th<tCount; th++) {
			AdvancementArg<T> *arg = new AdvancementArg<T>(cur, endSize, this, rnd);
			pthread_create(&threads[th], NULL, PopulatePool, (void*)arg);
		}
		next->reserve(endSize);
		PopulatePool(cur, next, rnd);
	
		for (uint th=0; th<tCount; th++) {
			void *rtn;
			pthread_join(threads[th], &rtn);
			AdvancementArg<T> *arg = (AdvancementArg<T> *)rtn;
			next->insert(next->end(), arg->newPop.begin(), arg->newPop.end());
			delete arg;
		}
		if (next->size() > endSize)
			next->resize(endSize);
		popSize = next->size();
		typename Population<T>::Type *t = cur;
		cur = next;
		next = t;
	}
	if (next != sourcePopulation) {
		Purge(next);
		delete next;
	}
	this->curPopulation = cur;
	return generationCount;
}

template <class T>
inline
void CPool<T>::WriteLdData(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh ) {
	assert(0);
	cerr<<"-----------------Error: WriteLdData() (cpool.h:555) has not been written. \n";

	return;
}

template<class T>
inline
void CPool<T>::BuildInitialPopulation(Utility::Random& rnd, int expressionCount) {
	/** 	
	LocusAssociationGrid *grid = NULL;
	This will be reintroduced if/when we incorporate the association matrix again
	if (gridSource.length() > 0) {
		grid = new LocusAssociationGrid();
		grid->Load(gridSource.c_str());
	} */
	if (sourcePopulation) {
		Purge(sourcePopulation);
		sourceAlleles.Clear();
	}
	else
		sourcePopulation = new vector<AlleleSource<T> *>();
	//sourcePopulation = new Population<T>::Type(expressionCount);
	int locusCount = loci->LocusCount();

	//If this fails, then we didn't initialize the loci vector
	assert(locusCount > 0);
	for (uint i=0; i<expressionCount; i++) {
		AlleleSourceSingle<T> *chrom = new AlleleSourceSingle<T>(loci, i);
		chrom->InitLoci(rnd, false);
		sourcePopulation->push_back(chrom);
		sourceAlleles.Add(chrom->GetSourceID(), chrom);
	}
		
/*	if (grid)
		delete grid; */
}

/*****/
template<class T>
inline
bool CPool<T>::BuildXOEvents(vector<size_t>& events, Utility::Random& rnd, float lambda) {
	if (lambda == 0.0)
		lambda = poissonLambda;
	size_t eventCount = PoissonEventCount(rnd, lambda);
	
	for (size_t i=0; i<eventCount; i++) {
		float loc = gdStart +  rnd(gdLength);
		Locus *locus = loci->At(loc);
		int idx = locus->GetID();
		vector<size_t>::iterator other = find(events.begin(), events.end(), idx);
		if (other != events.end())
			events.erase(other);
		else
			events.push_back(idx);
	}
	sort(events.begin(), events.end());
	return events.size() % 2;
}
/**
 * @brief when we are performing XO when gender is relevant, we adjust it according to (f=1.25N, m=0.75N) where N is the average lambda
 */
template<class T>
inline
bool CPool<T>::BuildXOEvents(vector<size_t>& events, Utility::Random& rnd, bool isMale) {
	float lambda;
	if (isMale)
		lambda = poissonLambda * 0.75;
	else
		lambda = poissonLambda * 1.25;
	return BuildXOEvents(events, lambda);
}	


template<class T>
inline
CPool<T>::CPool(int chromID, const char *project, const char *label) : curGeneration(0), isSuspended(false), id(chromID), project(project), label(label), sourcePopulation(NULL), curPopulation(NULL), loci(NULL), poissonLambda(0.0), gdStart(0.0), gdLength(0.0), remainingChroms(0), generationCount(0), sourceGeneration(0) { 
	pthread_mutex_init(&poolLock, NULL);
	pthread_mutex_init(&reportLock, NULL);
}

template<class T>
inline
CPool<T>::~CPool() {	
	if (curPopulation) {
		Purge(curPopulation);
		delete curPopulation;
	}

	if (sourcePopulation) {
		Purge(sourcePopulation);
		delete sourcePopulation;
	}
	pthread_mutex_destroy(&poolLock);
}


#ifdef CPPUNIT
class CPoolTest : public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( CPoolTest );
	CPPUNIT_TEST( TestInitialization );
	CPPUNIT_TEST_SUITE_END();
public:
	CPoolTest();
	~CPoolTest();

	void setUp();
	void tearDown();

	void TestInitialization();

};
#endif
}

#endif
