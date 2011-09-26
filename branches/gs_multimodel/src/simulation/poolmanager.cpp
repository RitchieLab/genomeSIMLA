//
// C++ Implementation: poolmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "poolmanager.h"
#include "basicsample.h"
#include "chromosome.h"
#include "utility/exception.h"

namespace Simulation {

using namespace StatusModel;
using namespace Utility::Exception;


uint PoolManager::simultaneousChrom			= 1;
uint PoolManager::threadsPerChrom			= 1;

bool PoolManager::generateOverviewLD		= true;
//uint PoolManager::maxThreadCount        	= 2;

pthread_mutex_t PoolManager::chromPoolLock 	= PTHREAD_MUTEX_INITIALIZER;	///<Lock down the pool whenever we are doing threaded stuff with it
pthread_mutex_t PoolManager::summaryLock	= PTHREAD_MUTEX_INITIALIZER;	///<Lock down the summary stream

//bool PoolManager::useDistanceMaps			= true;
//bool PoolManager::distanceMappingFunction	= NULL;


using namespace PopulationGrowth;

PoolManager::PoolManager() : currGeneration(0), poolIsPopulated(false), threadedPoolIdx(0){
#ifdef USE_XY
	poolxy = NULL;
	locusXY = NULL;
#endif //USE_XY
}

size_t PoolManager::GetPoolCount() {
#ifdef USE_XY
	if (poolxy)
		return pools.size() + 1;
#endif //USE_XY
	return pools.size();
}


PoolManager::~PoolManager() {
	size_t count = pools.size();
	for (uint i=0; i<count; i++) 
		delete pools[i];
#ifdef USE_XY
	if (poolxy)
		delete poolxy;

	if (locusXY)
		delete locusXY;
#endif //USE_XY
}


void PoolManager::AddPool(ChromPool *pool) {
	pools.push_back(pool);
}

string PoolManager::GenerateMarkerInfoFilename(uint pos, uint chr, const char *prj, uint startLoci, bool forDataset) {
	assert(chr < pools.size());
	stringstream ss;
	string dsPortion = "";
	if (forDataset)
		dsPortion = ".dataset";
#ifdef USE_XY
	if (chr == (uint)-1)
		ss<<prj<<"."<<pos<<"."<<poolxy->GetLabel()<<dsPortion<<"."<<startLoci<<".info";
	else
#endif //USE_XY
		ss<<prj<<"."<<pos<<"."<<pools[chr]->GetLabel()<<dsPortion<<"."<<startLoci<<".info";
	return ss.str();
}

string PoolManager::GenerateMarkerInfoFilename(uint gen, uint chr, const char *prj, bool forDataset) {
	assert(chr < pools.size());
	stringstream ss;
	string dsPortion = "";
	if (forDataset)
		dsPortion = ".dataset";
#ifdef USE_XY
	if (chr == (uint)-1)
		ss<<prj<<"."<<gen<<"."<<poolxy->GetLabel()<<dsPortion<<".info";
	else
#endif //USE_XY
		ss<<prj<<"."<<gen<<"."<<pools[chr]->GetLabel()<<dsPortion<<".info";
	return ss.str();
}



string PoolManager::GeneratePhasedFilename(uint gen, uint pos, const char *prj, bool forDataset) {
	assert(pos < pools.size());
	stringstream ss;
	string dsPortion = "";
	if (forDataset)
		dsPortion = ".dataset";
	ss<<prj<<"."<<gen<<"."<<pools[pos]->GetLabel()<<dsPortion<<".phased";
	return ss.str();
}

string PoolManager::GeneratePhasedFilename(uint gen, uint pos, const char *prj, int start, bool forDataset) {
	assert(pos < pools.size());
	stringstream ss;
	string dsPortion = "";
	if (forDataset)
		dsPortion = ".dataset";
	ss<<prj<<"."<<gen<<"."<<pools[pos]->GetLabel()<<dsPortion<<"."<<start<<".phased";
	return ss.str();
}



bool PoolManager::VerifyHaploview() {
	cout<<"Haploview verification is under construction\n";
	return true;
}

string PoolManager::GenerateSampleMarkerInfo(const char *prj, float minFreq) {
	stringstream filename;
	filename<<prj<<".dataset.marker";

	ofstream file(filename.str().c_str(), ios_base::out);

	/**
  	 * Iterate over each chromosome and write each locus (telling it what the freq threshold is)
	 */
	size_t count = pools.size();
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		ch->WriteMarkerInfo(file, minFreq);
	}
#ifdef USE_XY
	cerr<<"string PoolManager::GenerateSampleMarkerInfo(const char *prj, float minFreq)\n";
	cerr<<"This is just not a good idea! We are trying to squeeze locus information from non-XY and XY chromosomes into a single file! If we really need this to happen, we should consider changing the standard chromosome's output format to work reasonably well with the XY stuff\n";
	assert(false);

#endif //USE_XY

	return filename.str();
}



/**
 * @brief Generate phased output from each of the chromosomes as one file/chromosome
 */
string PoolManager::WritePhasedOverview(const char *project, bool forDataset) {
	size_t count = pools.size();
	float threshold = 0.0;
	string outputFilename = "";

	if (forDataset)
		threshold = Individual::minAlFreqThreshold;


	//Work through the various chromosome pools:
	// 1 Produce the overview image, if requested
	// 2 Work through the windows (increment by stride and show windows size number of snps 
	//   per window
	// 3 Produce the locus reports
	// * The phased output is always removed after we haploview generates the png
	// * The main overview pool is left behind if leavePhased is true
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		string filename = GeneratePhasedFilename(currGeneration, i, project, forDataset);

		ch->WriteBinaryPhased( filename.c_str());

		outputFilename = filename;
	}
#ifdef USE_XY
	if (poolxy)
		poolxy->Save(currGeneration);
#endif //USE_XY

	return outputFilename;
}


string PoolManager::WriteMarkerInfoOverview(const char *project, bool forDataset) {
	bool success = true;
	size_t count = pools.size();
	float threshold = 0.0;
	string outputFilename = "";

	if (forDataset)
		threshold = Individual::minAlFreqThreshold;


	//Work through the various chromosome pools:
	// 1 Produce the overview image, if requested
	// 2 Work through the windows (increment by stride and show windows size number of snps 
	//   per window
	// 3 Produce the locus reports
	// * The phased output is always removed after we haploview generates the png
	// * The main overview pool is left behind if leavePhased is true
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		string filename = GenerateMarkerInfoFilename(currGeneration, i, project, forDataset);
	
		ofstream file(filename.c_str(), ios_base::out);
		ch->WriteMarkerInfo( file, threshold );
	
		if (!file.is_open()){
			cout<<"Unable to write to chromosome data file: "<<filename<<". Verify that there is enough disk space and permissions are set properly. No LD maps were generated for this \n";
			success=false;
			break;
		}

		file.close();
		outputFilename = filename;
	}
#ifdef USE_XY
	if (poolxy) 
		OOOOPS("poolmanager.cpp:224");
#endif //USE_XY

	return outputFilename;
}



ChromPool *PoolManager::GetNextPool(uint &idx) {
	ChromPool *pool = NULL;
	POOLLOCK;

	if (threadedPoolIdx < pools.size())  {
		idx = threadedPoolIdx;
		pool = pools[threadedPoolIdx++];
	}
	POOLUNLOCK;
	return pool;
}


struct LdThread {
	bool doLoadFirst;
	PoolManager *mgr;
	ostream *summaryStream;
	ostream *locStream;
	LdThread(bool doLoadFirst, PoolManager *mgr, ostream *summaryStream, ostream *locStream) : doLoadFirst(doLoadFirst), mgr(mgr), summaryStream(summaryStream), locStream(locStream) { }
};
	
//This will be run as an independant thread (as well as the main thread)
//Each will pop one pool from the pools and work on them independantly until all
//are done
void *PoolManager::ProduceOverviewLD(void *arg) {
	uint idx = 0;
	stringstream summary;
	char filename[1024];

	LdThread thArg = *((LdThread *)arg);
	SUMMARYLOCK;
	//LocusReport locReport(thArg.mgr->locusReport);
	SUMMARYUNLOCK;

	uint currGen = thArg.mgr->GetCurrentGeneration();
	ChromPool *pool = thArg.mgr->GetNextPool(idx);
	
	if (pool) {
	
		while (pool) {
			sprintf(filename, "%s.%d.%s.ld", thArg.mgr->projectName.c_str(), currGen, pool->GetLabel().c_str());
			if (thArg.doLoadFirst) {
				string archiveFilename = thArg.mgr->GeneratePhasedFilename(currGen, idx, thArg.mgr->projectName.c_str(), false);
				pool->LoadArchive(archiveFilename.c_str());
				pool->WriteLdData(filename, summary, thArg.mgr->locusReport);
				SUMMARYLOCK;
				cout<<pool->GetLabel().c_str()<<" LD completed.\n";
				SUMMARYUNLOCK;
				pool->ClosePool(true);
			}
			else
				pool->WriteLdData(filename, summary, thArg.mgr->locusReport);
			pool = thArg.mgr->GetNextPool(idx);
		}

		
		SUMMARYLOCK;
		*(thArg.summaryStream)<<summary.str();
		thArg.mgr->locusReport.WriteHtmlReport( *(thArg.locStream) );
		//locReport.WriteHtmlReport( *(thArg.locStream) );

		SUMMARYUNLOCK;
	}
	return NULL;
}


string PoolManager::ProduceOverviewLD(const char *project, bool forDataset, bool loadFirst, ofstream &summary, const char *growthChart) {
	char locusReportFilename[1024];
	string sampleDetails;

	if (ChromPool::fastLD)
		sampleDetails = ".sampled";

	sprintf(locusReportFilename, "%s%s.%d.loc.html", project, sampleDetails.c_str(), GetCurrentGeneration());

	locusReport.Reset();

	summary<<"<P><CENTER><A HREF='"<<ExtractFilename(locusReportFilename)<<"'>Locus Selection Report</A>\n";
	summary<<"<P><CENTER><A HREF=\""<<ExtractFilename(growthChart)<<"\">Growth Chart</A>\n";

	ofstream lr(locusReportFilename, ios::out);
	lr<<"<HTML><HEADER><TITLE>Locus Report: Generation "<<GetCurrentGeneration();
	if (ChromPool::fastLD)
		lr<<"(Sampled)";
	lr<<"</TITLE></HEADER>\n";

	if (ChromPool::fastLD) 
		lr<<"<P><I>** This report is based on a sample drawn from a large population. Any Allele frequencies presented within the HTML document itself reflect the frequencies of the sample, not the entire population. In order to find out the frequency of alleles over the entire pool, please see the locus report at the top of the main report page.</I>\n<P>";
	

	uint maxThreadCount = simultaneousChrom - 1;

	pthread_t threads[maxThreadCount]; 

	//This has to be reset so we can start at the beginning
	threadedPoolIdx = 0;				

	LdThread thArgs(loadFirst, this, &summary, &lr);
	//Launch a thread for N-1 process where N is the max number of threads
	for (size_t i=0; i<maxThreadCount; i++) {
		pthread_create(&threads[i], NULL, ProduceOverviewLD, (void*)&thArgs);
	}
	//Do some work ourselves
	ProduceOverviewLD((void*)&thArgs);


#ifdef USE_XY
		//XY stuff is so different, we're cheating for now, and not trying to work it into the threading model
		if (poolxy) {
			char filename[32768];
			sprintf(filename, "%s.%d.%s.ld", projectName.c_str(), currGeneration, poolxy->GetLabel().c_str());
			poolxy->Refresh(currGeneration);
			poolxy->WriteLdData(filename, summary, locusReport, 0.0);
			cout<<"XY LD Completed\n";
			poolxy->Suspend();
		}
#endif //USE_XY

	//Let's catch up before we move on
	for (size_t i = 0; i<maxThreadCount; i++) {
		pthread_join(threads[i], NULL);
	}


	
	return locusReportFilename;
}



bool PoolManager::Dump(const char *prj) {
	bool success = true;
	uint count = pools.size();

	cout<<"Dumping the gene pool!\n";

	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		string filename = GeneratePhasedFilename(currGeneration, i, prj, false);
		ch->WriteBinaryPhased( filename.c_str() );

		ch->CalculateAlleleFrequencies(0);
		ch->SaveLoci(prj);
	}
#ifdef USE_XY
	if (poolxy) {
		poolxy->CalculateAlleleFrequencies(0);
		poolxy->Save(currGeneration);
	}
#endif //USE_XY
	return success;
}			

ChromPool *PoolManager::SeedChromosome(const char *seed, const char *loc, const char *label, const char *fmt, uint headerCount) {
	ChromPool *ch = new ChromPool(pools.size(), seed, loc, fmt, headerCount);
	ch->SetLabel( label );
	pools.push_back(ch);
	return ch;
}

ChromPool *PoolManager::AddChromosome(const char *locfilename, const char *gridFilename, const char *label) {
	ChromPool *ch = new ChromPool(pools.size(), locfilename, gridFilename);
	ch->SetLabel( label );
	pools.push_back(ch);
	return ch;
}
ChromPool *PoolManager::AddChromosome(const char *locFilename, const char *label) {
	ChromPool *ch = new ChromPool(pools.size(), locFilename);
	ch->SetLabel( label );
	pools.push_back(ch);
	return ch;
}

#ifdef USE_XY

void PoolManager::AddChromosomeXY(const char *locFilename, const char *label, Utility::Random& rnd) {
	if (poolxy) {
		cerr<<"Previous XY pool exists. Unwilling to load more than one XY\n";
		exit(1);
	}
	locusXY = new LocusManagerFileBased<LocusXY>(projectName.c_str(), -1);
	((LocusManagerFileBased<LocusXY>*)locusXY)->Load(locFilename);
	
	poolxy = new CPoolXY(99, projectName.c_str(), label);
	poolxy->Init(locusXY, rnd);
}


AlleleSource<LocusXY> *PoolManager::DrawX() {
	AlleleSource<LocusXY> *chrom = NULL;

	if (poolxy)
		chrom = poolxy->DrawX(Utility::Random::globalGenerator);
	return chrom;
}

AlleleSource<LocusXY> *PoolManager::DrawY() {
	AlleleSource<LocusXY> *chrom = NULL;

	if (poolxy)
		chrom = poolxy->DrawY(Utility::Random::globalGenerator);
	return chrom;

}
#endif //USE_XY
ChromPool *PoolManager::GetChromosome(uint idx) {
	ChromPool *pool = NULL;

	if (idx < pools.size()) 
		pool = pools[idx];
	return pool;
}

void PoolManager::DeleteChromosome(uint idx) {
	assert(idx < pools.size());
	ChromPool *ch = pools[idx];
	delete ch;

	vector<ChromPool *>::iterator itr = pools.begin();

	while (idx-- > 0)  
		itr++;
	
	pools.erase(itr);
}

ChromPool *PoolManager::AddChromosome(uint blockCount, ChromPool::BlockDefinition &defaultBlock, const char *label) {
	ChromPool *ch = new ChromPool(pools.size(), blockCount);
	ch->SetLabel( label );
	ch->DefineDefaultBlock( defaultBlock );
	pools.push_back(ch);
	return ch;
}

	
bool PoolManager::Close() {
	for (uint i=0; i<pools.size(); i++) {
		ChromPool *ch = pools[i];
		ch->ClosePool(false);
	}	
	poolIsPopulated=false;
	return true;
}

void PoolManager::InsertModel(DiseaseModel *model) {
	size_t poolCount = GetPoolCount();
	vector<uint> relLoci;
#ifdef USE_XY
	if (poolxy) {
		model->GetDiseaseLoci(-1, relLoci);
		poolxy->SetModelLoci(relLoci);
		poolCount--;
	}
#endif //USE_XY
	for (size_t i=0; i<poolCount; i++) {
		relLoci.clear();
		model->GetDiseaseLoci( i, relLoci);
		ChromPool *pool = pools[i];
		pool->ResetModelLoci();
		pool->SetModelLoci(&relLoci);
	}
}

		

bool PoolManager::PrepForSampling(const char *project, uint generation, DiseaseModel *model) {
	bool success = true;
	if (!poolIsPopulated) {
		size_t poolCount = GetPoolCount();
#ifdef USE_XY
		if (poolxy) 
			//I don't think we actually need to do anything to the new pools
			poolCount--;
#endif //USE_XY

		//BitSetType requiredChroms = model->AssociatedChromosomes(poolCount);

		vector<uint> relLoci;
		projectName = project;
		currGeneration = generation;
		for (uint i=0; success && i<poolCount; i++) {
			relLoci.clear();
			model->GetDiseaseLoci( i, relLoci);
			string filename = GeneratePhasedFilename(generation, i, project, false);
			ChromPool *pool = pools[i];
			pool->ResetModelLoci();
			pool->SetModelLoci(&relLoci);
			pool->ReadMinimalBinary( filename.c_str());
		}
	}
	return success;
}

bool PoolManager::ResolveGenotypes(vector<Individual*>& people, DiseaseModel *model) {
	bool success = true;
	if (!poolIsPopulated) {
		size_t poolCount = GetPoolCount();
#ifdef USE_XY 
		// The new pool object requires much less moderation, so this isn't necessary...lets just skip it
		if (poolxy)
			poolCount--;
#endif
		BitSetType chromsAlreadyInMem = model->AssociatedChromosomes(poolCount);
	
		for (uint i=0; success && i<poolCount; i++) {
			bool doLoadChrom = !chromsAlreadyInMem[i];
			doLoadChrom = true;
			ChromPool *ch = pools[i];
		
			if (doLoadChrom) {
				string filename = GeneratePhasedFilename(currGeneration, i, projectName.c_str(), false);
				ch->ReadBinaryPhased(filename.c_str(), true);
			}
			//	success = success && LoadPoolContents(projectName.c_str(), currGeneration, i, false);
			vector<Individual*>::iterator itr = people.begin();
			vector<Individual*>::iterator end = people.end();

			while (itr != end) {
				ch->ResolveGenotypes(*itr);
				itr++;
			}
			if (doLoadChrom) 
				pools[i]->ClosePool(true);
		}
	}
	return success;
}

bool PoolManager::ResolveGenotypes(const char *project, uint generation, Sample *samplePopulation, DiseaseModel *model) {
	bool success = true;
	if (!poolIsPopulated) {
		size_t poolCount = GetPoolCount();
#ifdef USE_XY 
		// The new pool object requires much less moderation, so this isn't necessary...lets just skip it
		if (poolxy)
			poolCount--;
#endif
		BitSetType chromsAlreadyInMem = model->AssociatedChromosomes(poolCount);
	
		for (uint i=0; success && i<poolCount; i++) {
			bool doLoadChrom = !chromsAlreadyInMem[i];
			doLoadChrom = true;
			if (doLoadChrom) 
				success = success && LoadPoolContents(project, generation, i, false);
		
	/*		cout<<"Pool ["<<i<<"]"<<"\n";
			pools[i]->SamplePhased(cout, 0, 25, 0, 25);
			cout<<"Sample of 25 items from pool "<<i<<"\n";
	*/	
			if (success)
				samplePopulation->ResolveGenotypes(pools[i], i, doLoadChrom);
			if (doLoadChrom) 
				pools[i]->ClosePool(true);
		}
	}
	return success;
}

void PoolManager::SetProjectName(const char *project) {
	projectName = project;
}

bool PoolManager::LoadLoci(const char *project, uint generation) {	
	bool success = true;
	projectName = project;
	uint poolCount = pools.size();

	for (uint i=0; success && i<poolCount; i++) {
		success = success && LoadLoci(project, generation, i);
	}
	

#ifdef USE_XY
		if (locusXY) {	
			((LocusManagerFileBased<LocusXY>*)locusXY)->SetPrefix(project);
			((LocusManagerFileBased<LocusXY>*)locusXY)->Refresh(generation);
			success = success && true;
		}
#endif //USE_XY
	if (success)
		currGeneration = generation;

	return success;
}

void PoolManager::SetPoolGeneration(uint generation) {
#ifdef USE_XY
	if (poolxy) 
		if (generation > 0)
			poolxy->Refresh(generation);
#endif //USE_XY
	uint poolCount = pools.size();
	locusMap.clear();
	for (uint i=0; i<poolCount; i++) 
		pools[i]->SetGeneration(generation);
	currGeneration = generation;
}

bool PoolManager::LoadLoci(const char *project, uint generation, uint chrID) {
	ChromPool *ch = pools[chrID];
	projectName = project;
	ch->SetGeneration(generation);
	bool success = true;
	try {
		if (ch->LoadLoci(project, locusMap, false) > 0)
			currGeneration = generation;
		else 
			success = false;
	} catch (Utility::Exception::FileNotFound& e) {
		success = false;
		cout<<e.GetErrorMessage()<<"\n";
	}
	return success;
}

bool PoolManager::LoadPoolContents(const char *project, uint generation, uint chrID, bool closeAfterLoading) {
	ChromPool *ch = pools[chrID];
	string filename = "";
	try {
		filename = GeneratePhasedFilename(generation, chrID, project, false);
	}
	catch (Utility::Exception::FileNotFound& e) {
		cout<<"Exception: "<<e.GetErrorMessage()<<"\n";
		return false;
	}
	ch->ReadBinaryPhased( filename.c_str(), true);
//	cout<<"Recalculating Allele Frequencies..";cout.flush();
//	ch->CalculateAlleleFrequencies();
//	cout<<"Done!\n";
	if (closeAfterLoading)
		ch->ClosePool(true);
	return true;
}

bool PoolManager::LoadPoolContents( const char *project, uint generation, bool closeAfterLoading ) {
	bool success = true;
	uint poolCount = pools.size();
	for (uint i=0; success && i<poolCount; i++) {
		success = success && LoadPoolContents(project, generation, i, closeAfterLoading);
	}
	if (closeAfterLoading)
		poolIsPopulated = false;
	if (success)
		currGeneration = generation;

	return success;
}		

bool PoolManager::ForceAlleleFrequency(uint chrID, uint locusID, float al1, float al2) {
	forcedFrequencies.push_back(ForcedAF(chrID, locusID, al1, al2));
	return true;
}


bool PoolManager::ForceAlleleFreqRange(uint chrID, uint locusID, float min, float max) {
	if (max < min)  
		return 0;

	float af1, af2;
	float diff = max-min;
	float value = min + Utility::Random::globalGenerator.drand() * diff;
	
	if (min > 0.5) {
		af2 = 1.0 - value;
		af1 = value;
	}
	else {
		af2 = 1.0 - value;
		af1 = value;
	}


	forcedFrequencies.push_back(ForcedAF(chrID, locusID, af1, af2));
	return true;
}


void PoolManager::DrawIndividual(Individual& ind, bool isXX) {
	int size = pools.size();
	ind.SetGender(isXX?2:1);
	for (int i=0; i<size; i++)  {
		pools[i]->DrawIndividual( &ind );
	}
#ifdef USE_XY
	if (poolxy)  {
		poolxy->DrawIndividual(ind, isXX, Utility::Random::globalGenerator);
	}
#endif	
}

void PoolManager::DrawIndividual(Individual& ind) {
	bool isXX = Utility::Random::globalGenerator.drand() < 0.5;
	DrawIndividual(ind, isXX);
}
bool PoolManager::GetAlleleFrequency(uint chrID, uint locId, float &al1, float &al2) {
	bool success = false;
	
#ifdef USE_XY
	if (chrID == -1)	
		if (poolxy)
			success = poolxy->GetAlleleFrequencies(locId, al1, al2);
		else
			throw Utility::Exception::General("One or more disease loci were assigned to the X/Y chromosome, but none has been configured");
	else 
#endif
	if (chrID < pools.size() ) 
		success = pools[chrID]->GetAlleleFrequencies(locId, al1, al2);
	return success;
}

void PoolManager::CalculateAlleleFrequencies() {
	int size = pools.size();
	for (int i=0; i<size; i++)  
		pools[i]->CalculateAlleleFrequencies(0);
}

struct GenThread {
	bool doLoadFirst;
	PoolManager *mgr;
	uint moreGenerations;
	GrowthRate *f;

	GenThread(bool doLoadFirst, PoolManager *mgr, uint moreGenerations, GrowthRate *f) : doLoadFirst(doLoadFirst), mgr(mgr), moreGenerations(moreGenerations), f(f) { }
};


void *PoolManager::AdvanceGenerations(void *args) {
	GenThread thArg = *((GenThread*)args);
	
	ChromPool::continueRunning = true;
	uint idx = 0;
	ChromPool *ch = thArg.mgr->GetNextPool(idx);
	string filename;

	while (ch && ChromPool::continueRunning) {
		int curGen = ch->GetGenerationCount();
		if (thArg.doLoadFirst)  {
			//uint resultGeneration = thArg.mgr->currGeneration+thArg.moreGenerations;
			filename = thArg.mgr->GeneratePhasedFilename(thArg.mgr->currGeneration, idx, thArg.mgr->projectName.c_str(), false);
			
			try {
				ch->LoadArchive( filename.c_str() );
			} catch (Exception::FileIO &e) {
				if (thArg.mgr->currGeneration == 0) {
					cout<<"Starting at generation 0, but unable to load. Rebuilding the loci\n";
				}
			}
		}
		uint resultGeneration = ch->AdvanceGenerations( thArg.moreGenerations, thArg.f, threadsPerChrom);

		filename = thArg.mgr->GeneratePhasedFilename(resultGeneration, idx, thArg.mgr->projectName.c_str(), false);

		ch->WriteBinaryPhased( filename.c_str());
		ch->SaveLoci(thArg.mgr->projectName.c_str());
		if (thArg.doLoadFirst)
			ch->ClosePool(true);

		SUMMARYLOCK;
		cout<<"\n"<<ch->GetLabel().c_str()<<" advanced "<<ch->GetGenerationCount() - curGen<<" (pop. "<<ch->GetExpressionCount()<<").\n";
		SUMMARYUNLOCK;

		ch = thArg.mgr->GetNextPool(idx);
	}

	return NULL;
}

bool PoolManager::AdvanceGenerations(uint moreGenerations, GrowthRate *f, bool loadFirst) {
	bool success = true;
	
	threadedPoolIdx = 0;

	GenThread thArgs(loadFirst, this, moreGenerations, f);

	uint maxThreadCount = simultaneousChrom - 1;
	pthread_t threads[maxThreadCount]; 

	for (size_t i=0; i<maxThreadCount; i++) 
		pthread_create(&threads[i], NULL, AdvanceGenerations, (void*)&thArgs);

	AdvanceGenerations((void*)&thArgs);
	
	for (size_t i=0; i<maxThreadCount; i++) 
		pthread_join(threads[i], NULL);

#ifdef USE_XY
	//Cheating by not incorporating the XY stuff into the threading model
	if (poolxy) {
		int curGen = poolxy->GetCurrentGeneration();
		poolxy->Wake();
		curGen = poolxy->AdvanceGenerations(moreGenerations, f, Utility::Random::globalGenerator, simultaneousChrom);
		poolxy->Bake();
		poolxy->CalculateAlleleFrequencies(poolxy->GetPopulationSize());
		poolxy->Save(curGen);
		((LocusManagerFileBased<LocusXY>*)locusXY)->Save(curGen, poolxy->GetLabel().c_str());
		poolxy->Suspend();
		currGeneration = poolxy->GetCurrentGeneration();
	}
	else
#endif //USE_XY
	//Changed since we might exit at target generation
		currGeneration = pools[0]->GetGenerationCount();
	//currGeneration+=moreGenerations;
	return success;
}

size_t PoolManager::GetExpressionCount() {
#ifdef USE_XY
	if (poolxy) 
		return poolxy->GetPopulationSize();
#endif //USE_XY
	assert(pools.size() > 0);
	return pools[0]->GetExpressionCount();
}

void PoolManager::CreateInitialPopulation(uint populationCount, bool doDump,  bool closeAfterCreating) {
#ifdef USE_XY
	if (poolxy) {
		poolxy->BuildInitialPopulation(Utility::Random::globalGenerator, populationCount);
		poolxy->Suspend();
	}
#endif //USE_XY
	uint count = pools.size();
	doDump = doDump || closeAfterCreating;			///<We want to make sure we dump if we are about to close
	for (uint i=0; i<count; i++) {
		ChromPool *pool = pools[i];
		pool->BuildPool(populationCount);
	
		string filename;
		//Binary Phased Output
		if (doDump || closeAfterCreating) {
			filename = GeneratePhasedFilename(currGeneration, i, projectName.c_str(), false);
			pool->WriteBinaryPhased( filename.c_str() );

			if (closeAfterCreating)
				pool->ClosePool(true);
		}
		
		pool->SaveLoci(projectName.c_str());
	}
	poolIsPopulated = true;
}




void PoolManager::GenerateReport(ostream &os, uint headerWidth) {
	size_t count = pools.size();
	os<<setw(headerWidth)<<"Chromosome Pool Count: "<<count<<endl;
	os<<setw(headerWidth)<<"Chromosome Intialization Method:";
	if (ChromPool::UseAdamEve) 
		os<<"2 Founder Chromosomes\n";
	else if (ChromPool::UseEden)
		os<<"2 Founders Crossed to Produce Initial Population\n";
	else 
		os<<"Randomized\n";

	for (uint i=0; i<count; i++) {
		os<<setw(headerWidth)<<pools[i]->GetLabel()<<"\n";
		os<<setw(headerWidth + 5)<<"Loci Count: "<<pools[i]->GetLociCount()<<"\n";
		if (pools[i]->GetExpressionCount() > 0)
			os<<setw(headerWidth + 5)<<"Expression Count: "<<pools[i]->GetExpressionCount()<<"\n";
		os<<setw(headerWidth + 5)<<"XO Lambda: "<<pools[i]->XOLambda()<<"\n";
	}
}

bool PoolManager::GetLocusReference(const char *label, Locus& locus) {
	//Protect against misreading end of lines
	if (strlen(label) == 0)
		return false;
	bool success = locusMap.find(label) != locusMap.end();
	if (success) 
		locus = locusMap[label];
	else {
		LocusMap::iterator itr = locusMap.begin();
		LocusMap::iterator end = locusMap.end();

		cout<<" Can't find "<<label<<" in our map ("<<locusMap.size()<<"). Here is what we have\n";
		while (itr != end) {
			cout<<"Locus "<<itr->first<<" : "<<itr->second.GetChromID()<<":"<<itr->second.GetLabel()<<"\n";
			itr++;
		}
	}
	return success;
}

bool PoolManager::InitializeLoci(uint startGen, const char *prj, bool doLoad) {
	projectName = prj;


	//Create our own loci if we are starting at the beginning
	uint startID = 0;
	if (!doLoad && startGen == 0) {
		//cout<<"Clearing locusMap ! PoolManager::InitializeLoci() 680\n";
		locusMap.clear();
		vector<ChromPool*>::iterator itr = pools.begin();
		vector<ChromPool*>::iterator end = pools.end();
		while (itr != end) {
			uint nextStart = (*itr)->InitLoci(startID, locusMap);
			if (nextStart == startID) {
				vector<ChromPool*>::iterator prev = itr;
				itr++;
				pools.erase(prev);
				delete *prev;
			} else {
				startID = nextStart;
				itr++;
			}
		}
		for (uint i=0; i<forcedFrequencies.size(); i++) {
			ForcedAF f=forcedFrequencies[i];
#ifdef USE_XY
			if (f.chrID == (uint)-2) 
				locusXY->ForceAlleleFrequency(f.locus, f.af1, f.af2);
			else 
#endif	
			if (f.chrID >= pools.size() || !pools[f.chrID]->ForceAlleleFrequency(f.locus, f.af1, f.af2, locusMap)) 
				cout<<"\nAn error was encountered trying to force allele frequencies: "<<f.chrID<<" "<<f.locus<<" "<<f.af1<<" "<<f.af2<<"\n";
		}
		poolIsPopulated = pools.size() > 0;
	}
	else {
		poolIsPopulated = LoadLoci(prj, startGen);

	}
#ifdef USE_XY
		if (locusXY) {	
			locusXY->BuildLocusMap(locusMap);
		}
#endif //USE_XY

	return poolIsPopulated;
}

bool PoolManager::InitializePools(uint startGen, uint poolSize, const char *prj, 
			bool doLoad, bool doDump, bool closeAfterCreating) {
//	size_t count = pools.size();
	bool success = false;
	projectName = prj;
#ifdef USE_XY
		if (locusXY)
			((LocusManagerFileBased<LocusXY>*)locusXY)->SetPrefix(prj);
#endif	
	//Create our own loci if we are starting at the beginning
//	uint startID = 0;
	if (!doLoad && startGen == 0) {
//cout<<"Clearing locusMap InitializePools("<<startGen<<","<<poolSize<<","<<prj<<","<<doLoad<<","<<doDump<<","<<closeAfterCreating<<") 703\n";
//locusMap.clear();
		CreateInitialPopulation( poolSize, doDump, closeAfterCreating );
		poolIsPopulated = !closeAfterCreating;
		success = true;
	
	} 
	else {
		//SetPoolGeneration(startGen);
//		if (!closeAfterCreating)
		success = LoadPoolContents(prj, startGen, closeAfterCreating);
//		else
#ifdef USE_XY
		if (locusXY)
			((LocusManagerFileBased<LocusXY>*)locusXY)->SetPrefix(prj);
		if (poolxy) {
			poolxy->SetProject(prj);
			poolxy->Open(locusXY, Random::globalGenerator);
			poolxy->Refresh(startGen);
		}
#endif			
	}

	return success;

}

void PoolManager::SetLocusSelectors(LocusSelectionArray& locSel) { 
	LocusSelectionArray::iterator itr = locSel.begin();
	LocusSelectionArray::iterator end = locSel.end();

	while (itr != end) {
		locusReport.AddLocusSelection( *itr++ );
	}
}



}
