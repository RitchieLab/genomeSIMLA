//Chrompol.cpp

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

#include <sstream>
#include <iomanip>
#include <math.h>
#include "chrompool.h"
#include "utility/random.h"
#include "ldplotter.h"
#include "utility/strings.h"
#include "genetics/allelicdecoder.h"
#include "ldpngcomponent.h"
#include "utility/exception.h"

using namespace Genetics;

namespace Simulation {

using namespace Visualization;

float 	ChromPool::defFre1 					= 0.5;
float 	ChromPool::defFre2 					= 0.5;
float 	ChromPool::errRate 					= 0.01;
bool 	ChromPool::randomizeAlleleFreq 		= true;
float 	ChromPool::maxLdConsideration 		= 500000;
size_t 	ChromPool::numberOfBlocksToReport 	= 15;
//DistMappingFn *ChromPool::mappingFn 		= new KosambiMapping();
size_t 	ChromPool::LdSpread					= 250;
bool 	ChromPool::writeDPrimePlots 		= true;
bool 	ChromPool::writeRSquaredPlots 		= true;
uint 	ChromPool::highBufferSize 			= 50;
uint 	ChromPool::detailsBufferSize 		= 5;
string  ChromPool::cssFilename				= "http://chgr.mc.vanderbilt.edu/genomeSIMLA/genomesim.css";
bool 	ChromPool::closePoolBetweenAdvances = false;
bool 	ChromPool::UseOriginalCrossing		= false;
bool 	ChromPool::UseAltLD					= true;


bool 	ChromPool::writeLdTextReport		= false;		
uint	ChromPool::reportWidth				= 1250;
bool	ChromPool::fastLD					= true;
uint 	ChromPool::maxLDIndividuals			= 1500;
uint	ChromPool::plotScanSize				= 1000;
uint 	ChromPool::minOffspring				= 1;
uint 	ChromPool::maxOffspring				= 1;
uint 	ChromPool::targetPop 				= 0;
bool	ChromPool::continueRunning			= true;
bool 	ChromPool::UseAdamEve				= false;
bool 	ChromPool::UseEden					= false;
int		ChromPool::FounderCount				= 25;
double	ChromPool::FounderDistortion		= 0.25;
int		ChromPool::MaxRepeatCount			= 5;
double	ChromPool::ParentDistortion			= 0.0;
double	ChromPool::ChildDistortion			= 0.0;
bool	ChromPool::WritePhasedPools			= false;
size_t 	ChromPool::PhasedPoolWriteSize		= 1500;


pthread_mutex_t ChromPool::poolLock 		= PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t ChromPool::reportLock		= PTHREAD_MUTEX_INITIALIZER;

using namespace std;
using namespace Utility;

void ChromPool::SetOffspringLimits(uint min, uint max) {
	minOffspring = min;
	maxOffspring = max;
	if (max < min) {
		throw Exception::General("Offsprint limits for max and min don't make sense.");
	}
	cout<<"Setting generational advancement to produce between "<<min<<" and "<<max<<" offspring per coupling. This does not strictly follow hardy-weinburg pricipals\n";
}
//
// Default constructor
//
	//Sleep requirements
	bool neverSleep;					///<Marked true if part of a disease model. 
	bool isSleeping;					
	string archiveFilename;				///<Make sure we know where to find the pool when we wake up

ChromPool::ChromPool(uint chromID, uint blockCount) : reloadLocusFile(true),
	type(BlockBased), chromID(chromID), generationCount(0), lociCount(0), blockCount(blockCount), locSource(""), gridSource(""), seedSource(""), summaryReport(NULL), isMinimal(false),seedFmt(""), headerCount(0) {
	frequenciesDetermined = false;
	positionsDefined = false;
}
ChromPool::ChromPool(uint chromID, const char *locFilename) : reloadLocusFile(true),type(FileBased), chromID(chromID), generationCount(0), 
	lociCount(0), blockCount(0), locSource(locFilename), gridSource(""), seedSource(""), summaryReport(NULL), isMinimal(false), seedFmt(""), headerCount(0) {
	frequenciesDetermined = false;
	positionsDefined = false;
}

ChromPool::ChromPool(uint chromID, const char *locFilename, const char *gridFilename) : reloadLocusFile(true),type(FileBased), chromID(chromID), generationCount(0), 
	lociCount(0), blockCount(0), locSource(locFilename), gridSource(gridFilename), seedSource(""), summaryReport(NULL), isMinimal(false), seedFmt(""), headerCount(0) {
	frequenciesDetermined = false;
	positionsDefined = false;
}
ChromPool::ChromPool(uint chromID, const char *seedFilename, const char *locFilename, const char *fmt, uint headerCount) : reloadLocusFile(true),
	type(FileBased), chromID(chromID), generationCount(0), lociCount(0), blockCount(0), 
	locSource(locFilename), gridSource(""), seedSource(seedFilename), summaryReport(NULL), isMinimal(false), 
	seedFmt(fmt), headerCount(headerCount) {

	frequenciesDetermined = false;
	positionsDefined = false;
}

bool ChromPool::UseMapInfoFile() {
	return positionsDefined;
}


/**
 * @NOTE The following functions, BuildPool_Eden and BuildPool_AdamEve were
 * used for experimentation, and aren't considered to be suitable for 
 * generating valid chromosome pools. 

void ChromPool::BuildPool_Eden(uint expressionCount) {
	Chromosome adam(&loci, poissonLambda, &recombIndexLookup);
	Chromosome eve(&loci, poissonLambda, &recombIndexLookup);
	adam.InitLoci();
	cout<<"** Eden Based Pool Initialization\n";
	cout<<"Founder Count:            "<<FounderCount<<"\n";
	cout<<"Founder Distortion Level: "<<FounderDistortion<<"\n";
	cout<<"Max Repeats:	             "<<MaxRepeatCount<<"\n";
	cout<<"Parent Distortion:        "<<ParentDistortion<<"\n";
	cout<<"Child Distortion:         "<<ChildDistortion<<"\n";
	
	vector<Chromosome> initPop;
		
	for (int i=0; i<FounderCount; i++) {
		adam.Distort(FounderDistortion);
		int repeatCount = rnd.lrand(1, MaxRepeatCount-1);

		for (int n=0; n<repeatCount; n++) {
			initPop.push_back(adam);
		}
		eve = adam;
		eve.Invert();
		//eve.Distort(FounderDistortion);	
		repeatCount = rnd.lrand(1, MaxRepeatCount-1);

		for (int n=0; n<repeatCount; n++)	{
			initPop.push_back(eve);
		}
	}
	int founderCount = initPop.size() - 1;
	cout<<"Total Founder Count:       "<<founderCount+1<<"\n";
	//If this fails, then we didn't initialize the loci vector
	assert(loci.size() > 0);
	size_t eventCount;
	if (seedSource.length() == 0) {
		for (uint i=0; i<expressionCount; i++) {
			Chromosome ch(&loci, poissonLambda, &recombIndexLookup);
			int idx1=rnd.lrand(0, founderCount);
			int idx2=rnd.lrand(0, founderCount);
			Chromosome mom = initPop[idx1];
			Chromosome dad = initPop[idx2];
			mom.Distort(ParentDistortion);
			dad.Distort(ParentDistortion);
			ch = mom.CrossImmediate(dad, eventCount);
			ch.Distort(ChildDistortion);
			pool.push_back(ch);
		}	
	}
	else 	
		LoadSeedData(expressionCount);
	
}
 */

/*
/// Initializes the the pool with expressionCount/2 copies of Adam and expressionCount/2 copies of eve
void ChromPool::BuildPool_AdamEve(uint expressionCount) {
	Chromosome founder(&loci, poissonLambda, &recombIndexLookup);
	founder.InitLoci();
	cout<<"** Adam/Eve Based Pool Initialization\n";
	cout<<"Founder Count:            "<<FounderCount<<"\n";
	cout<<"Founder Distortion Level: "<<FounderDistortion<<"\n";
	cout<<"Max Repeats:	             "<<MaxRepeatCount<<"\n";
	cout<<"Child Distortion:         "<<ChildDistortion<<"\n";
	
	vector<Chromosome> initPop;
		
	for (int i=0; i<FounderCount; i++) {
		Chromosome newFounder(&loci, poissonLambda, &recombIndexLookup);
		newFounder = founder;
		//newFounder = founder;
		newFounder.Distort(FounderDistortion);
		int repeatCount = rnd.lrand(1, MaxRepeatCount-1);
		for (int n=0; n<repeatCount; n++)
			initPop.push_back(newFounder);

		newFounder.Invert();
		newFounder.Distort(FounderDistortion);
	
		repeatCount = rnd.lrand(1, MaxRepeatCount-1);
		for (int n=0; n<repeatCount; n++)
			initPop.push_back(newFounder);
	}
	int founderCount = initPop.size() - 1;
	cout<<"Total Founder Count:       "<<founderCount<<"\n";
	//If this fails, then we didn't initialize the loci vector
	assert(loci.size() > 0);

	if (seedSource.length() == 0) {
		for (uint i=0; i<expressionCount; i++) {
			Chromosome ch(&loci, poissonLambda, &recombIndexLookup);
			int idx = rnd.lrand(0,founderCount);
			ch = initPop[idx];
			ch.Distort(ChildDistortion);
			pool.push_back(ch);
		}
	}
	else 	
		LoadSeedData(expressionCount);
}
*/
	
double ChromPool::XOLambda() {
	return poissonLambda;
}
void ChromPool::BuildPool(uint expressionCount) {
	LocusAssociationGrid *grid = NULL;
	if (gridSource.length() > 0) {
		grid = new LocusAssociationGrid();
		grid->Load(gridSource.c_str());
	}
/*	if (UseAdamEve) {
		BuildPool_AdamEve(expressionCount);	
		return;
	}
	else if (UseEden) {
		BuildPool_Eden(expressionCount);
		return;
	}
*/
	//If this fails, then we didn't initialize the loci vector
	assert(loci.size() > 0);
	if (seedSource.length() == 0) {
		for (uint i=0; i<expressionCount; i++) {
			Chromosome ch(&loci, poissonLambda, &recombIndexLookup);
			ch.InitLoci(grid);
			pool.push_back(ch);
		}
	}
	else 	
		LoadSeedData(expressionCount);
		
	if (grid)
		delete grid;

}

/**
 * @brief Add a potential block of loci the chromosome
 * @param minCount The minimum number of snps 
 * @param maxCount The the maximum number of snps
 * @param minMap The minimum recombination Fraction between one locus and and the next one
 * @param maxMap The max recombination fraction between one locus and the next one
 * @param blockStrength Adjusts the block away from 0.5 chance. The larger the value, the weaker the "block" will be
 */
uint ChromPool::DefineBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax, float frequency) {
	assert(minSnpCount <= maxSnpCount);
	float min=0.000001;				///<Weird roundoff, this allows us to accept 0.000001
	if (blckMin < min || blckMax < min || snpMin < min || snpMax < min) {
		throw Utility::Exception::General("Block definitions should be defined in Centi-Morgans, and should be 1/100000 or greater. ");
	}

	blockPrototypes.push_back(BlockDefinition(minSnpCount, maxSnpCount, blckMin, blckMax, snpMin, snpMax, frequency, blockPrototypes.size() + 1));
	return blockPrototypes.size();
}


void ChromPool::ClearBlocks(){ 
	blockPrototypes.clear();
}


ChromPool::BlockDefinition &ChromPool::GetDefaultBlock() {
	return defaultBlock;
}

ChromPool::BlockDefinition &ChromPool::GetBlock(size_t idx) {
	size_t size = blockPrototypes.size();
	assert(idx < size);

	return blockPrototypes[idx];
}


size_t ChromPool::GetBlockCount() {
	return blockPrototypes.size();
}

void ChromPool::DefineDefaultBlock(ChromPool::BlockDefinition& block) {
	defaultBlock = block;
}

void ChromPool::DefineDefaultBlock(uint minSnpCount, uint maxSnpCount, float blckMin, float blckMax, float snpMin, float snpMax) {
	assert(minSnpCount <= maxSnpCount);
	defaultBlock = ChromPool::BlockDefinition(minSnpCount, maxSnpCount, blckMin, blckMax, snpMin, snpMax, 100.0, 0);
}

ChromPool::BlockDefinition &ChromPool::DrawBlockDefinition() {
	float draw = Utility::Random::globalGenerator.drand();
	bool blockFound = false;
	uint blockCount = blockPrototypes.size();
	BlockDefinition &theBlock = defaultBlock;
	float currValue=0.0;
	for (uint i=0; i<blockCount && !blockFound; i++) {
		currValue += blockPrototypes[i].frequency;
		if (draw <= currValue) {
			theBlock = blockPrototypes[i];
			blockFound = true;
		}
	}
	return theBlock;		
}

double ChromPool::GetMapDistance(uint locus) {
	double rec = 0.0;
	if (locus<lociCount)
		rec = loci[locus].MapDistance();
	return rec;
}




//This format can't support missing data at the moment!
void ChromPool::WriteBinaryPhased(const char *filename) {
	ofstream file(filename, ios::out|ios::trunc|ios::binary);
	WriteBinaryPhased(file);
	file.close();


	if (WritePhasedPools) {
		stringstream pedFilename;
		pedFilename<<filename<<".ped";
		file.open(pedFilename.str().c_str(), ios::out|ios::trunc);
		SaveAsPhased(file, 0.0, 0, PhasedPoolWriteSize);
		file.close();
	}
}


void ChromPool::SamplePool(size_t expCount, size_t genotypeCount) {
	
	size_t poolSize = pool.size();
	cout<<"Pool Sample: "<<chromID<<" has "<<poolSize<<" expressions\n";
	if (poolSize < expCount)
		expCount = poolSize;

	for (size_t i=0; i<expCount; i++) {
		cout<<i<<"\t: ";
		pool[i].ShowGenotypes(genotypeCount, ' ');
		cout<<"\n";
	}
	
	if (poolSize > expCount) {
		cerr<<"...\n";
		for (size_t i=poolSize - expCount; i<poolSize; i++) {
			cerr<<i<<"\t";	
			pool[i].ShowGenotypes(genotypeCount, ' ');
			cerr<<"\n";
		}
	}
}


void ChromPool::WriteBinaryPhased(ofstream &file) {
	size_t expressionCount = GetExpressionCount();
	file.write((char*)&chromID, 4);
	file.write((char*)&generationCount, 4);
	file.write((char*)&expressionCount, 4);
	assert(pool.size() > 0);
	int lc = pool[0].LociCount();
	file.write((char*)&lc, 4);

	for (size_t chID=0; chID<expressionCount; chID++) {
		pool[chID].WriteBinary(file);
	}
}

void ChromPool::ReadMinimalBinary(const char *filename) {


	//Make sure we don't waste time loading a pool for nothing
	if (modelLoci.size() == 0) 
		return;

	minPool.clear();
	poolDetails.underlyingLoci = &loci;
	poolDetails.poissonLambda = poissonLambda;
	poolDetails.recombIndexLookup = &recombIndexLookup;

	ifstream file(filename, ios::in|ios::binary);
	if (file.is_open()) {
		uint expressionCount = GetExpressionCount();
		file.read((char*)&chromID, 4);
		file.read((char*)&generationCount, 4);
		file.read((char*)&expressionCount, 4);
		uint lociCount;
		file.read((char*)&lociCount, 4);
		
		minPool.clear();
		pool.clear();
		for (size_t chID=0; chID<expressionCount; chID++) {
			Chromosome c(&loci, poissonLambda, &recombIndexLookup);
			c.ReadBinary(file, lociCount, &modelLoci, false);
//			MinChrom minchrom(&modelLoci, &poolDetails);
//			minchrom = c;
			pool.push_back(c);
			//pool.push_back(c);
			//pool[chID] = c;
		}

	}
	else {
		throw Utility::Exception::FileNotFound(filename);
	}
	isMinimal = true;
}

void ChromPool::ReadBinaryPhased(const char *filename, bool forceLoad) {

	if (!forceLoad && pool.size() > 0) 
		return;
	ifstream file(filename, ios::in|ios::binary);
	if (file.is_open()) 
		ReadBinaryPhased(file);
	else {
		throw Utility::Exception::FileNotFound(filename);
	}
		
}


void ChromPool::ResetModelLoci() {
	modelLoci.clear();
}

void ChromPool::SetModelLoci(vector<uint> *ml) {
	modelLoci = *ml;
}

void ChromPool::AddModelLocus(Locus &l) {
	modelLoci.push_back(l.GetID());
}

void ChromPool::ReadBinaryPhased(ifstream &file) {
//	int maxPrints = 25;	

	size_t expressionCount = GetExpressionCount();
	file.read((char*)&chromID, 4);
	file.read((char*)&generationCount, 4);
	file.read((char*)&expressionCount, 4);
	size_t lociCount = 0;
	file.read((char*)&lociCount, 4);

	pool.clear();
	pool.resize(expressionCount);
	for (size_t chID=0; chID<expressionCount; chID++) {
		Chromosome c(&loci, poissonLambda, &recombIndexLookup);
		c.ReadBinary(file, lociCount, &modelLoci, true);
	
		/*if (maxPrints-- > 0) {
			cout<<chID<<" --> ";
			c.(10, ' ');
			cout<<"\n";
		}*/
		pool[chID] = c;
	}
}

uint ChromPool::WriteMarkerInfo(std::ostream& os, float minAlFreq /*=0*/) {
	return WriteMarkerInfo(os, 0, lociCount, minAlFreq);
}

void ChromPool::WriteMendelFormat(ostream& os) {
	LocusArray::iterator itr = loci.begin();
	LocusArray::iterator end = loci.end();

	while (itr != end)
		itr++->WriteMendelFormat(os);
	
}

uint ChromPool::WriteMarkerInfo(std::ostream& os, uint first, 
		uint count, float minAlFreq) {
	if (first == lociCount)
		return 0;

	uint last=0;
	
	if (count == 0)
		last = lociCount - first;
	else {
		last = first+count;
		if (last > lociCount)
			last = lociCount;
	}

	for(uint i=first; i<last; i++) {
		Locus &l = loci[i];
		if (l.GetAlleleFreq(1) >= minAlFreq && l.GetAlleleFreq(2) >= minAlFreq) 
			l.WriteMarkerInfo(os);
	}

	return last;
}

//We have to sort the loci by their position, since the map position is dependant on sorted loci
struct locPosLT {
	bool operator()(const Locus& left, const Locus& right) {
		return left.GetLocation() < right.GetLocation();
	}
};
bool ChromPool::SortLoci() {
	//Let's keep a copy of the unsorted items for now, so we can determine if they were out of order
	LocusArray sortedLoci = loci;
	
	sort(sortedLoci.begin(), sortedLoci.end(), locPosLT());

	bool changesRequired = false;
	
	for (uint i=0; i<loci.size(); i++) {
		Locus &sl = sortedLoci[i];
		Locus &l = loci[i];
		changesRequired = changesRequired || sl.GetID() != l.GetID();
		sl.SetID(i);
	
		if (changesRequired) {
			cout<<GetLabel()<<" "<<sl.GetID()<<" "<<sl.GetLabel()<<" "<<sl.GetLocation()<<" :\t: "<<l.GetID()<<" "<<l.GetLabel()<<"  "<<l.GetLocation()<<"\n";
		}
	}

	loci.clear();
	loci = sortedLoci;
	
	return changesRequired;
}

//bool ChromPool::ResolveGenotypes(
void ChromPool::ClearLoci() {
	loci.clear();
}

void ChromPool::ValidateBlocks() {
	bool validBlockConfigurations = defaultBlock.minBlckMap > 1e-15 && defaultBlock.maxBlckMap > 1e-15 && defaultBlock.minSnpMap > 1e-15 && defaultBlock.maxSnpMap > 1e-15;
	if (validBlockConfigurations) {
		size_t blockCount = blockPrototypes.size();
		while (blockCount--) {	
			BlockDefinition &block = blockPrototypes[blockCount];
			if (block.minBlckMap < 1e-15 || block.maxBlckMap < 1e-15 || block.minSnpMap < 1e-15 || block.maxSnpMap < 1e-15) {
				char s[1028];
				sprintf(s, "Block %d has invalid recombination bounds.", block.id);
				throw Utility::Exception::General(s);
			}		
		}
	}
	else	{
		string s = "The default block has invalid recombination bounds. ";
		throw Utility::Exception::General(s.c_str());
	}
}

uint ChromPool::InitLoci(uint starterID, LocusMap &locusMap) {
	int markerID = 0;	//starterID;
	recombIndexLookup.Clear();


	if (locSource.length() == 0) {

		ValidateBlocks();
		int lPos = 1;					///<The physical location of the locus

		for (uint i=0; i<blockCount; i++){ 
			BlockDefinition &block = DrawBlockDefinition();
			uint snpCount = block.GetSnpCount();
			float map = block.GetBlockMapDistance();
			if (map < 0.000001)
				map = 0.000001;

			stringstream blckID;
			blckID<<"Blck"<<i<<"x"<<block.id;

			for (uint snpID=0; snpID<snpCount; snpID++) {
				lPos += (int)(map * 1000000);

				Locus l(map, chromID, markerID, lPos);
				l.SetLabel(chromID, ++markerID + starterID);
				if (snpID == 0) 
					l.SetDescription(blckID.str().c_str());

				//Determine how we assign allele frequencies
				if (randomizeAlleleFreq)
					l.RandomizeFreq(defFre1, defFre2);
				else
					l.AssignFreq(defFre1, defFre2);
				//Add the loci to the collection
				loci.push_back(l);
//				cout<<"* Locus: "<<l.GetLabel()<<" "<<snpID<<" : "<<l.GetMinAlleleFreq()<<"\n";
				//Adjust the recombination value for the next locus
				//recomb=minRecomb+(rDiff * Utility::Random::globalGenerator.drand());
				map = block.GetSnpMapDistance();
//				cout<<"\tlocusMap["<<l.GetLabel()<<"]="<<l.GetDescription()<<"\n";
				locusMap[l.GetLabel()] = l;
			}
		}
		lociCount = loci.size();
		if (SortLoci())
			cout<<GetLabel()<<" was reordered. Please be sure that any locus references uses an updated locus report for disease assignment or SNP id (RS Number) to ensure agreement on SNP selection for models. \n";
	
		float mapUnits = 0.0;
		for (size_t i=0; i<lociCount; i++) {
			Locus &l = loci[i];
			mapUnits+=l.MapDistance();
			recombIndexLookup.Insert(mapUnits, i);
			l.MapPosition(mapUnits);
		}
	
		poissonLambda=(mapUnits - loci[0].MapDistance()) * 0.01;
	}
	else {
		ifstream lf(locSource.c_str(), ios_base::in);
		if (lf.is_open()) {
			if (LoadLoci(lf, locusMap, true) == 0)
				return starterID;
			else
				markerID = loci[loci.size()-1].GetLocation();
		} 
		else {
			cout<<"Unable to open locus file '"<<locSource<<"'.\n";
			lociCount = 0;
			throw Utility::Exception::FileNotFound(locSource.c_str());
		}
	}

	return markerID + starterID;
}

void ChromPool::ClosePool(bool keepMinimalGenotypes) {
	minPool.clear();


	poolDetails.underlyingLoci = &loci;
	poolDetails.poissonLambda = poissonLambda;
	poolDetails.recombIndexLookup = &recombIndexLookup;
	
	if (keepMinimalGenotypes) {
		PoolType::iterator itr = pool.begin();
		PoolType::iterator end = pool.end();
		while (itr != end) {
			MinChrom minchrom(&modelLoci, &poolDetails);
			minchrom = *itr;
			minPool.push_back(minchrom);
			itr++;
		}	
		isMinimal = true;
	} else {
	}
	pool.clear();

}


string ChromPool::GenerateLocusFilename(const char *prj) {
	stringstream ss;
	ss<<prj<<"."<<generationCount<<"."<<GetLabel()<<".loc";
	cout<<ss.str()<<"\n";
	return ss.str();
}


void ChromPool::SaveLoci (const char *project, float minThreshold) {
	string filename = GenerateLocusFilename(project);

	ofstream file(filename.c_str(), ios_base::out );
	if (file.is_open()) {
		SaveLoci( file, minThreshold );
		file.close();
	}
	else {	
		cout<<"FileIOError: "<<filename<<"\n";
		throw Utility::Exception::FileNotWritable(filename.c_str());
	}	
}


uint ChromPool::LoadLoci(const char *project, LocusMap& lmap, bool forceLoad) {
	string filename = locSource;
	if (generationCount) 
		filename = GenerateLocusFilename(project);
	bool success = false;

	ifstream file(filename.c_str(), ios_base::in);
	if (file.is_open()) {
		success = LoadLoci(file, lmap, forceLoad);
		file.close();
		reloadLocusFile = false;
	}
	else {
		throw Utility::Exception::FileNotFound(filename.c_str());
	}
	return success;

}

uint ChromPool::LoadLoci(istream& data, LocusMap& locusMap, bool forceLoad /* = false */) {
	frequenciesDetermined = true;

	if (!reloadLocusFile && !forceLoad && loci.size() > 0)
		return loci.size();

	lociCount = 0;
	
	char line[4096];
	data.getline(line, 4096);			///<The chromosome id
	data.getline(line, 4096);			///<The count of loci
	data.getline(line, 4096);			///<The header information

	//string id;

	//We need to make sure we are starting with an empty locus pool
	loci.clear();
	while (!data.eof()) {
		Locus loc(chromID, lociCount);
		//data>>id;
		if (!data.eof()) {
			data>>loc;

			if (loc.Valid()) {
				loci.push_back(loc);
				lociCount++;			///<Only want to increment if the locus is to be kept
			}

			//For debugging
			//cout<<id<<" "<<loc<<" \n";
			
		}
	}
	lociCount = loci.size();
	if (SortLoci())
		cout<<GetLabel()<<" was reordered. Please be sure that any locus references uses an updated locus report for disease assignment or SNP id (RS Number) to ensure agreement on SNP selection for models. \n";

	//Well, let's check the last locus for it's position to determine if there are valid positions 
	if (lociCount > 0)
		positionsDefined = loci[lociCount - 1].GetLocation() > 0;

	float mapUnits = 0.0;			///<Total genetic distance covered
	for (size_t i=0; i<lociCount; i++) {
		Locus &l = loci[i];
 		mapUnits+= l.MapDistance();
		l.MapPosition(mapUnits);

		//Make sure we keep no duplicates
		recombIndexLookup.Set(mapUnits, i);
		
		int attempt = 1;
		string curLabel = l.GetLabel();
		while (locusMap.find(l.GetLabel()) != locusMap.end()) {
			stringstream ss;
			ss<<curLabel<<"-"<<attempt++;
			l.SetLabel(ss.str().c_str());
		}
		if (attempt > 1)
			cout<<"Duplicate RS numbers were encountered: "<<curLabel<<" was renamed to "<<l.GetLabel()<<".\n";
		
		//cout<<"\t+ locusMap["<<l.GetLabel()<<"] = "<<l.GetDescription()<<" ("<<locusMap.size()<<")\n";
		locusMap[l.GetLabel()] = l;
	}

	poissonLambda=(mapUnits - loci[0].MapDistance()) * 0.01;

	//We don't count the first distance, since it specifies the distance from itself and a non-existant previous SNP
//	for (size_t i=1; i<loci.size(); i++ )
//		poissonLambda+=loci[i].MapDistance();
	reloadLocusFile = false;
	return lociCount;



}



void ChromPool::LoadPhased(istream& data, uint currGen, uint headerCount, bool forceLoad /* = false */) {

	if (!forceLoad && pool.size() > 0)
		return;

	//We certainly don't want to double our pool size!
	pool.clear();

	uint expectedLineLength = ((lociCount * 2) + 2) * 2;
	char *line = new char[expectedLineLength];
	data.getline(line, expectedLineLength);
	data.seekg(0, ios::beg);

	generationCount = currGen;

	uint columnCount = CountColumns((const char*)line);
	if (columnCount != lociCount + 2) {
		throw ("wrong number of SNPs!");
		//cout<<"The current gene pool doesn't match the configuration. Please check correct the settings or rerun the initial population once again.\n\tLoci Count in File: "<<((float)(columnCount-6)/2.0)<<" - Configuration file specifies: "<<lociCount<<"\n";
	}

	while (!data.eof()) {
		string var;
		
		//Ignore the first 6 columns
		for (uint i=0; i<2; i++) 
			data>>var;

		if (!data.eof() ) {
			string gbg;						///<Just a place holder for garbage
			//Get the two strands of DNA 
			Chromosome c(&loci, poissonLambda, &recombIndexLookup);
			for (uint i=0; i<headerCount; i++) 
				data>>gbg;

			for (uint i=0; i<lociCount; i++) {
				int strand1;
				data>>strand1;
				c.SetValue( i, strand1 );
			}
			pool.push_back(c);
		}
	}
	
	lociCount = loci.size();
	delete[] line;
	cout<<"Previous gene pool has been loaded. Configuration parameters might not explicitly describe the current set of individuals. \n";
}

void ChromPool::LoadSeedData(uint expressionCount) {
	ifstream seedData(seedSource.c_str());
	if (!seedData.is_open()) {
		cout<<"Unable to open seed datafile: "<<seedSource<<"\n";
		abort();
	}
	AllelicDecoder decoder("?");

	//For now, we'll assume a line is no more than 104K characters wide
	char line[1048576];
	seedData.getline(line, 1048576);

	while (!seedData.eof()) {
		seedData.getline(line, 1048576);
		decoder.SeedEncodings(line, 1, 4);
	}
	decoder.Verify(2);

	//For some reason, I'm not having any luck with seekg and tell
//	while (pool.size()) {
		seedData.close();
		seedData.open(seedSource.c_str());
		seedData.getline(line, 1048576);
		vector<int> individual;
		while (!seedData.eof()) {
			seedData.getline(line, 1048576);
			if (decoder.ParseLine(line, individual, 1, 4, 0) ) {
				//translate this into a bitset	
	
				Chromosome chrom(&loci, poissonLambda, &recombIndexLookup);
				chrom.InitLoci(individual);
	//			BitSetType missing(count, false);
				pool.push_back(chrom);
				
			}
		}
//	}	
	cout<<"Seed Data Contained: "<<pool.size()<<" individuals.\n";
/*	if (pool.size() + ((float)pool.size() * 1.05) < expressionCount) {
		cout<<"The seed data in file "<<seedSource<<" doesn't have enough individuals to meet the initial population. Please adjust your growth parameters\n";
		abort();
	}
	*/
	CalculateAlleleFrequencies();
}

bool ChromPool::LoadArchive(const char *filename) {
	ReadBinaryPhased(filename, isMinimal);

	if (pool.size() == 0) {
		cout<<"Pool archive, "<<filename<<", contained no pool data. Please make sure there is enough diskspace to write the pool data to file\n";
		abort();
	}
	isMinimal = false;
	return true;
}
/*
uint ChromPool::AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f, const char *filename) {
	ReadBinaryPhased(filename, false);
	if (pool.size() == 0)  {
		cout<<"Pool archive, "<<filename<<", contained no pool data. Please make sure there is enough diskspace to write the pool data to file\n";
		abort();
	}
	
	uint poolSize = AdvanceGenerations(genCount, f);
	ClosePool(true);
	
	return poolSize;
}
*/
struct AdvancementArg {
	ChromPool::PoolType &sourcePop;
	ChromPool::PoolType newPop;
	uint count;
	ChromPool *pool;

	AdvancementArg(ChromPool::PoolType &sourcePop, uint count, ChromPool *pool) : sourcePop(sourcePop), count(count), pool(pool) { 
		newPop.reserve(count);
	}
};


uint ChromPool::AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f, uint threadCount) {
	threadCount--;

	uint goalPop = targetPop;
	if (goalPop == 0)
		goalPop = (uint)-1;

	PoolType newPool;
	
	PoolType *cur = &pool;
	PoolType *next = &newPool;

	uint popSize = pool.size();
	frequenciesDetermined = false;
	for(uint currGen=0; currGen < genCount && popSize < goalPop && continueRunning; currGen++){   
		next->clear();
	
		uint endSize = (*f)(++generationCount);

		if (endSize > goalPop)
			endSize = goalPop;

		servedPools = endSize;
		cout<<".";cout.flush();
		//uint portions = endSize / (threadCount + 1) + (5 * (threadCount + 1));
		//uint served = 0;
		pthread_t threads[threadCount];
		for (uint th=0; th<threadCount; th++) {
			AdvancementArg *arg = new AdvancementArg(*cur,endSize, this);
			pthread_create(&threads[th], NULL, PopulatePool, (void*)arg);
			//served+=portions;
		}
		next->reserve(endSize);	
		
		PopulatePool(*cur, *next);


		for (uint th=0; th<threadCount; th++) {
			void *rtn;
			pthread_join(threads[th], &rtn);
			AdvancementArg *arg = (AdvancementArg *)rtn;
			next->insert(next->end(), arg->newPop.begin(), arg->newPop.end());
			delete arg;
		}
		if (next->size() > endSize)
			next->resize(endSize);
		popSize = next->size();
		PoolType *t = cur;
		cur = next;
		next = t;
  	}
	this->pool = *cur;

	usedIndividuals.clear();
	return generationCount;
}


void ChromPool::WriteConfiguration(ofstream &oFile) {
	if (type == BlockBased) {
		oFile<<"ADD_CHROMOSOME "<<blockCount<<" "<<label<<"\n";
		uint bc = blockPrototypes.size();
		for (size_t i=0; i<bc; i++) {
			blockPrototypes[i].WriteConfiguration(oFile);
		}
	}
	else if (type == FileBased) {
		oFile<<"LOAD_CHROMOSOME \""<<locSource<<"\" "<<label<<"\n";
	}
	else {
		cout<<"Unknown chromosome type: "<<type<<"\n";
		abort();
	}
}

void *ChromPool::PopulatePool(void *args) {
	AdvancementArg *thArg = (AdvancementArg*)args;
	thArg->pool->PopulatePool(thArg->sourcePop, thArg->newPop);
	return (void*)args;
};



void ChromPool::PopulatePool(PoolType &source, PoolType &dest) {
	int origPop = source.size();

	size_t eventsObserved;						///<Number of X-Over events observed
//	int count = 0;
//	int totalEvents = 0;
	while (ContinueGrowing()) {
		int firstChromIndex = int(Random::globalGenerator(origPop));
		int secondChromIndex = int(Random::globalGenerator(origPop));

		assert(firstChromIndex < origPop);
		assert(secondChromIndex < origPop);
		//LOCKPOOL;

		//By default, this will be 1
//		int numOffspring = rnd.lrand((int)minOffspring, (int)maxOffspring);
//		while (numOffspring-- > 0 && doContinue) {
Chromosome c = source[firstChromIndex].CrossImmediate(source[secondChromIndex], eventsObserved);
//count++;
//totalEvents+=eventsObserved;
dest.push_back(c);
//			dest.push_back(source[firstChromIndex].CrossImmediate(source[secondChromIndex], eventsObserved));
//			doContinue = true;
//		}
		//UNLOCKPOOL;
	}
//cerr<<"pool completed. Average Events: "<<totalEvents/count<<"\n";

}


uint ChromPool::AdvanceGenerations(uint genCount, PopulationGrowth::GrowthRate *f) {
	int firstChromIndex, secondChromIndex;
	frequenciesDetermined = false;
	
	uint popSize = pool.size();
	cout<<"*";cout.flush();
	for(uint currGen=0; currGen < genCount && popSize < targetPop; currGen++){   
		PoolType lastpool = pool;
	 	int numChroms = lastpool.size();

		cout<<"Poolsize: "<<numChroms;
		pool.clear();
		uint endSize = (*f)(++generationCount);

		size_t eventsObserved;						///<Number of X-Over events observed

		for(uint currChrom=0; currChrom < endSize; currChrom++){
			
			firstChromIndex = int(Random::globalGenerator(numChroms));
			secondChromIndex = int(Random::globalGenerator(numChroms));
			if (UseOriginalCrossing)
				pool.push_back(lastpool[firstChromIndex].Cross(lastpool[secondChromIndex], eventsObserved));
			else {
				//cout<<"Event Count: "<<eventsObserved<<"\n";
				pool.push_back(lastpool[firstChromIndex].CrossImmediate(lastpool[secondChromIndex], eventsObserved));
			}
    	}
		popSize = pool.size();
		cout<<" - "<<popSize<<"\n";
  	}
	usedIndividuals.clear();
	return generationCount;
}




bool ChromPool::ForceAlleleFrequency(uint locusID, float al1, float al2, LocusMap &locusMap) {
	bool success = locusID < lociCount;

	if (success) {
		Locus &loc = loci[locusID];
		loc.AssignFreq(al1, al2);

		locusMap[loc.GetLabel()] = loc;
	}
	else {
		cout<<"Unable to set frequency for "<<GetLabel()<<":"<<locusID + 1<<"\n";
		cout<<"There is no locus "<<locusID<<" (Total loci "<<lociCount + 1<<")\n";
	}

	return success;
}



/**
 * 0 based
 */
bool ChromPool::GetAlleleFrequencies(uint locus, float &af1, float &af2) {
	if (locus<lociCount) {
		af1 = loci[locus].Freq1();
		af2 = loci[locus].Freq2();
	}
	
	return locus<lociCount;
}

size_t ChromPool::CalculateAlleleFrequencies(uint maxIndividualCount) {

	uint fixedAlleles = 0;
	uint indCount = pool.size();

	if (indCount == 0)
		return 0;
	if (maxIndividualCount == 0)
		maxIndividualCount = indCount;	size_t lociCount = loci.size();

	//locix2 array for each allele
	uint alleleCount = 2 * lociCount;
	//MMmmmm, MS won't bother to comply with C99, so I gotta go back to mid 90s 
	uint *freq = new uint[alleleCount];
	//uint freq[lociCount][2];

	for (uint i=0; i<alleleCount; i++) 
		freq[i]=0;

//Just testing something
	if (maxIndividualCount > 0 && indCount > maxIndividualCount) {
		indCount = maxIndividualCount;
	}

	for (uint ind = 0; ind<indCount; ind++) {
		Chromosome &curChrom = pool[ind];

		for(uint i=0; i<lociCount; i++) {
			freq[i*2+curChrom.At(loci[i].GetID())]++;
		}
	}

	uint idx = 0;
	for (uint i=0;i<alleleCount; ) {
		if (freq[i] == 0)
			fixedAlleles++;
		double al1 = (double)(freq[i++])/(double)(indCount);
		if (freq[i] == 0)
			fixedAlleles++;
		double al2 = (double)(freq[i++])/(double)(indCount);;
		loci[idx++].AssignFreq(al1, al2);
	}


	if (fixedAlleles == lociCount) {
		for (uint i=0; i<indCount; i++) {
			cout<<"Loci Count: "<<pool[i].LociCount()<<"\t";
			pool[i].ShowGenotypes(60, '>');
			cout<<"\n";
		}
		cout<<"We have a problem. There are no variants in our pool at all. Every one of our "<<indCount<<" expressions is the same\n";
		abort();
	}

	delete[] freq;
	return fixedAlleles;
}


void ChromPool::WriteStyleSheetDetails(ostream &os) {
/*	os<<"<STYLE type=\"text/css\">\n";
	os	<<"body { color:#000000; background-color: #ffffff; link:#918744; vlink:#84BA57; alink:#7381ba }\n"
		<<"h1 { border-width: 0; background-color:#eeeeee; color: #444C56; text-align: center }\n"
		<<"table { width: 95%; align: center; frame: border; rules: all }\n"
		<<"th { color: #eeeeee; background-color: #444c56 }\n"
		<<"th.invert { color:#444c56; background-color: #ffffff }\n"
		<<"</STYLE>\n";
*/
	os<<"<LINK REL='stylesheet' TYPE='text/css' HREF=\""<<ChromPool::cssFilename<<"\">\n";
}

void ChromPool::GenerateReport( ostream &os, uint headerWidth) {
	os<<setw(headerWidth)<<"Close Pools When not in Use: "; closePoolBetweenAdvances?os<<"Yes\n":os<<"No\n";
	os<<setw(headerWidth)<<"Draw DPrime: "; (writeDPrimePlots?os<<"Yes\n":os<<"No\n");
	os<<setw(headerWidth)<<"Draw RSquared: "; writeRSquaredPlots?os<<"Yes\n":os<<"No\n";
	os<<setw(headerWidth)<<"Write LD Report (text): "; writeLdTextReport?os<<"Yes\n":os<<"No\n";
	os<<setw(headerWidth)<<"Max Distance for Pairwise LD: "<<LdPlotter::maxSnpDistance<<"\n";
	os<<setw(headerWidth)<<"Num. Blocks Reported per Chromosome: "<<numberOfBlocksToReport<<"\n";
	os<<setw(headerWidth)<<"Perform Sampled LD: "; fastLD?os<<"Yes\n":os<<"No\n";
	os<<setw(headerWidth)<<"Max. Snps per Row: "<<ImageParameters::maxSnpsPerRow<<"\n";
	os<<setw(headerWidth)<<"LD Buffer Size: "<<LdSpread<<"\n";
	os<<setw(headerWidth)<<"Font: "<<ImageParameters::font<<"\n";
	os<<setw(headerWidth)<<"CSS Filename: "<<cssFilename<<"\n";
}

void ChromPool::InitSummaryReport(const char *filename) {
	stringstream ss;
	ss<<filename<<"-SUMMARY_REPORT.csv";
	summaryReport = new ofstream(ss.str().c_str(), ios::out|ios::trunc);
	*summaryReport<<"Generation\tPopulation Count\tFixed Alleles\tBlock Count\tMax Block Size\n";
}

void ChromPool::WriteLdData(const char *basefilename, ostream &summary, Visualization::LocusReport &locReport, float thresh ){
	char reportFilename[2048];
	stringstream summaryTemp;			///<Just used to hold onto reporting info until the end
	string filename = basefilename;
	if (fastLD)
		filename += "-sampled";

	sprintf(reportFilename, "%s.html", filename.c_str());
	ofstream blockReport(reportFilename, ios::out);
	blockReport<<"<HTML><HEADER><TITLE>Block Report for Generation "<<generationCount<<"</TITLE></HEADER>\n";
	blockReport<<"<BODY>\n";
	WriteStyleSheetDetails(blockReport);

	//If we are supposed to just do quickie LD for scanning a growth curve....
	uint indThresh = 0;
	if (fastLD)
		indThresh = maxLDIndividuals;

	uint fixedAlleles = CalculateAlleleFrequencies(indThresh);

//Just testing something
	vector<Locus> testLoci;

	if (fastLD) {
		uint start = 0;
		uint end = loci.size() - start;
	
		
		if (fastLD && (end > start + plotScanSize)) 
			end = start+ plotScanSize;
		
	
	
		for (uint i=start; i<end; i++) 
			testLoci.push_back(loci[i]);
	}
	else
		testLoci = loci;
//		Visualization::LdPlotter ld(thresh, pool, loci);

	Visualization::LdPlotter ld(thresh, pool, testLoci);

	if (summaryReport == NULL)
		InitSummaryReport(filename.c_str());

	if (UseAltLD)
		ld.CalculateLD2(filename.c_str());
	else
		ld.CalculateLD();
	if (!continueRunning) 
		return;

	//vector<HaplotypeBlock *> *haplotypes = ld.Init();
	BlockList *haplotypes = ld.Init();

	LOCKREPORT;
	locReport.Init(testLoci, haplotypes, ExtractFilename(reportFilename).c_str(), GetLabel().c_str());
	UNLOCKREPORT;

	//For now, let's use the LdTextReport
	uint count = haplotypes->size();

	sort(haplotypes->begin(), haplotypes->end(), BlockSizeEval());
	cout<<".";cout.flush();
	uint detailedReportCount = numberOfBlocksToReport;
	if (detailedReportCount == 0 || count > numberOfBlocksToReport)
		detailedReportCount = numberOfBlocksToReport;

	char ldPlotName[1024];

	sprintf(ldPlotName, "%s", GetLabel().c_str());
	ld.SetLabel(GetLabel().c_str());
	blockReport<<"<DIV class='image-frame'>\n\n<CENTER>\n<H1>Block Report For "<<GetLabel();
	if (fastLD) 
		blockReport<<" - Sample Based";
	blockReport<<"</H1>\n";
	blockReport<<"\t<TABLE><CAPTION>Chromosome Details</CAPTION><TR><TH>Threshold</TH><TH>Pool Size</TH><TH>Total Number<BR>of Snps</TH><TH>Fixed Alleles</TH><TH>Block Count</TH></TR>\n";
	blockReport<<"\t\t<TR class='invert'> <TD>"<<thresh<<"</TD> <TD>"<<pool.size()<<"</TD> <TD>"<<loci.size()<<"</TD><TD>"<<fixedAlleles<<"</TD> <TD>"<<count<<"</TD></TR>\n\t</TABLE>\n\n<DIV class='spacer'></DIV>\n";

	summaryTemp<<"\t\t<TR><TD>\n\t\t<TABLE><CAPTION>"<<GetLabel();
	if (fastLD) 
		summaryTemp<<" (Sample based LD) ";
	summaryTemp<<"</CAPTION>\n";

	summaryTemp<<"\t\t\t<TR               ><TD><B>SNP Count</B></TD><TD>"<<loci.size()<<"</TD></TR>\n";
	summaryTemp<<"\t\t\t<TR class='invert'><TD><B>Fixed Alleles</B></TD><TD>"<<fixedAlleles<<"</TD></TR>\n";
	summaryTemp<<"\t\t\t<TR               ><TD><B>Block Count</B></TD><TD>"<<count<<"</TD></TR>\n\t\t</TABLE>\n";
	summaryTemp<<"\t\t</TD>\n";

	uint imgHeight, imgWidth;

	if (writeDPrimePlots) {
		cout<<".";cout.flush();
		string imgFilename = ld.WriteLdDPrime(filename.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
		char *picFilename = strrchr((char*)imgFilename.c_str(), '/');
		if (picFilename)
			picFilename++;
		else
			picFilename = (char*)imgFilename.c_str();
	
		

		if (imgWidth > reportWidth)
			imgWidth = reportWidth;
	
		blockReport<<"\t<DIV><img src='"<<Utility::ExtractFilename(picFilename)<<"' width=95% onClick=\"window.open('"<<Utility::ExtractFilename(picFilename)<<"');\"></img></DIV>\n";

		summaryTemp<<"\t\t<TD><DIV><img src='"<<ExtractFilename(picFilename)<<"' width=90%  onClick=\"window.open('"<<Utility::ExtractFilename(reportFilename)<<"');\"></img></DIV></TD></TR>\n";
	}

	if (writeRSquaredPlots) {
		cout<<".";cout.flush();
		string imgFilename = ld.WriteLdRSquared(filename.c_str(), NULL, imgHeight, imgWidth, 0, "sum");
		const char *picFilename = strrchr(imgFilename.c_str(), '/') ;
		if (picFilename)
			picFilename++;
		else
			picFilename = imgFilename.c_str();
	
		if (imgWidth > reportWidth)
			imgWidth = reportWidth;
	
		blockReport<<"\t<DIV><img src='"<<Utility::ExtractFilename(picFilename)<<"' width=95% onClick=\"window.open('"<<Utility::ExtractFilename(picFilename)<<"');\"></img><DIV>\n";

		if (!writeDPrimePlots)
			summaryTemp<<"\t\t<TD><DIV><img src='"<<picFilename<<"' width=90%  onClick=\"window.open('"<<Utility::ExtractFilename(reportFilename)<<"');\"></img></DIV></TD></TR>\n";
	}

	if (writeLdTextReport) {
		cout<<".";cout.flush();
		string txtReportFilename = Utility::ExtractFilename(ld.WriteLdReport(filename.c_str(), NULL).c_str());
		blockReport<<"\t<P><A HREF='"<<txtReportFilename<<"' onClick=\"window.open('"<<txtReportFilename<<"');\">LD Values</A>\n";

	}

	LOCKREPORT;
	summary<<summaryTemp.str();
	UNLOCKREPORT;

	string locusFilename = Utility::ExtractBaseFilename(basefilename) + ".loc";
	blockReport<<"\t<P><A HREF='"<<locusFilename<<"' onClick==\"window.open('"<<locusFilename<<"');\">Locus Report</A>\n";
	blockReport<<"\t<P>"<<locusFilename<<"\n";

	

	*summaryReport<<generationCount<<"\t"<<pool.size()<<"\t"<<fixedAlleles<<"\t"<<count<<"\t";
	summaryReport->flush();
	if (count>0) {
		BlockListNode *block = (*haplotypes)[0];
		*summaryReport<<block->BlockCount()<<"\t";
		cout<<".";cout.flush();	
	}
	*summaryReport<<"\n";

	//We'll write the detailed reports here first. Then dump them to the file later
	stringstream details;
	string anchor;
	if (count > 0)
		(*haplotypes)[0]->HtmlSummaryBegin(blockReport);
	for (uint i=0;i<count && continueRunning; i++) {
		cout<<".";cout.flush();
		BlockListNode *block = (*haplotypes)[i];
		//HaplotypeBlock *block = (*haplotypes)[i];
		if (i<detailedReportCount) {
			//block->Report(cout);
			//block->DetailedReport(cout);
			sprintf(ldPlotName, "%s.%d", filename.c_str(), i+1);
			string dprimeFN, rsquaredFN, dprimeDet, rsquaredDet;
			if (writeDPrimePlots) {
				dprimeFN = ld.WriteLdDPrime(ldPlotName, block, imgHeight, imgWidth, highBufferSize, "-ovr");
				dprimeDet = ld.WriteLdDPrime(ldPlotName, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			}
			if (writeRSquaredPlots) {
				rsquaredFN = ld.WriteLdRSquared(ldPlotName, block, imgHeight, imgWidth, highBufferSize, "-ovr");
				rsquaredDet = ld.WriteLdRSquared(ldPlotName, block, imgHeight, imgWidth, detailsBufferSize, "-det");
			}
			
			anchor = block->HtmlDetailed(details, imgWidth, dprimeFN.c_str(), rsquaredFN.c_str(), dprimeDet.c_str(), rsquaredDet.c_str());
				
			block->HtmlSummary( blockReport, anchor.c_str() );
			//ld.WriteLdReport( ldPlotName, block, highBufferSize );	
		}
		else {
			block->HtmlSummary( blockReport, NULL );
			//block->Report(cout);
		}
	}
	if (count> 0)
		(*haplotypes)[0]->HtmlSummaryEnd( blockReport );
	blockReport<<details.str();
	blockReport<<"</BODY></HTML>\n";
	cout<<"\n";
}





void ChromPool::SaveLoci(ostream& os, float minThreshold) {
	CalculateAlleleFrequencies();
    ios_base::fmtflags old_settings = os.flags();
	os.setf(ios::fixed, ios::floatfield);

	int width=16;

	os << "Locus Log for "<<GetLabel()<<"\n";
	os << lociCount<<" Loci\n";

	if (lociCount > 0) {
		if (loci[0].GetMinAlleleFreq() >= minThreshold)
		//os << setw(width)<<"Locus ID" << " ";
			loci[0].WriteHeader(os, width);
	}
	for(uint i=0; i<lociCount; i++){
			//os << setw(width) << i+1 << " ";
		if (loci[i].GetMinAlleleFreq() >= minThreshold)
			loci[i].WriteFormatted(os, width);
	}
	os.setf(old_settings);

}


void ChromPool::SamplePhased(std::ostream& os, uint firstExpression, uint expCount, uint firstSnp, uint snpCount) {
	PoolType::iterator i;
	PoolType::iterator end=pool.end();

	uint lastExpression=firstExpression+expCount;
	if (lastExpression == 0 || lastExpression > pool.size())
		lastExpression = pool.size();

	uint lastSnp = firstSnp + snpCount;
	if (lastSnp == 0 || lastSnp > pool[0].LociCount())
		lastSnp = pool[0].LociCount();

	uint id=firstSnp;
	for (uint i=firstExpression; i<lastExpression; i++){ 
		os<<++id<<" 1 ";
		pool[i++].WritePedFormat(os, 0.0, firstSnp, lastSnp);

		os<<++id<<" 1 ";
		pool[i++].WritePedFormat(os, 0.0, firstSnp, lastSnp);
	}
}

//Let's dump the gene pool to a stream. We'll use haploview format
uint ChromPool::SaveAsPhased(std::ostream& os, float minThresh, uint first, uint count) {
	if (first == pool.size())
		return 0;

	PoolType::iterator i;
	PoolType::iterator end=pool.end();

	uint last=0;
	if (count == 0)
		last = pool[0].LociCount() - first;
	else {
		last = first+count;
		if (last > pool[0].LociCount())
			last = pool[0].LociCount();
	}


	uint id=first;
	for (i=pool.begin(); i!=end && i+1 != end; ) {
		os<<++id<<" 1 ";
		(*i++).WritePedFormat(os, minThresh, first, last);

		os<<id<<" 1 ";
		(*i++).WritePedFormat(os, minThresh, first, last);
	}

	return last;
}

void ChromPool::SaveAsPhased(std::ostream& os, float minThresh) {

	PoolType::iterator i;
	PoolType::iterator end=pool.end();

	uint id=0;

	for (i=pool.begin(); i!=end && i+1 != end; ) {
		//i is pedigree. 1 is individual ID. All the rest of metadata is 0
		os<<++id<<" 1 ";
		(*i++).WritePedFormat(os, minThresh, 0, 0);

		os<<id<<" 1 ";
		(*i++).WritePedFormat(os, minThresh, 0, 0);
	}
}

/*float ChromPool::GetErrorRate(uint locus) {
	assert(locus<lociCount);

	return loci[locus].ErrorRate();
}
*/

void ChromPool::ReserveIndividual(Individual *ind) {
	string id = ind->GetPoolID(chromID);
	//usedIndividuals[id]=ind->GetID();			//or something like this
	usedIndividuals[id]=1;
}



void ChromPool::ApplyChromosomeIDs(Individual *ind) {
	uint ChromCount = pool.size();
	if (ChromCount == 0)
		ChromCount = minPool.size();

	uint chromID1;
	uint chromID2;

	char id[128];

	//Make sure we don't redraw any individual
	IndLookupType::iterator pos;
	IndLookupType::iterator end = usedIndividuals.end();

	uint attempts = 0;

	do {
		chromID1 = int(Random::globalGenerator.drand() * ChromCount);
		chromID2 = int(Random::globalGenerator.drand() * ChromCount);
		sprintf(id, "%dx%d", chromID1, chromID2);
		pos = usedIndividuals.find(id);
		if (attempts++>5000) {
			cout<<"We are having trouble getting individual "<<ind->GetID()<<" from chromosome: "<<chromID<<". \nThe list of used individuals is "<<usedIndividuals.size()<<" long\n";
			cout<<"Go back to your data and make sure that your pool is big enough to support the number of datasets and individuals that are being requested. With your current settings, the pool is: "<<ChromCount<<"\n";
			abort();
		}
	} while (pos != end);

	ind->SetPoolID(chromID, id);
	ind->SetChromIDs(chromID, chromID1, chromID2);


}

Chromosome ChromPool::GetExpression(uint chrIDX) {
	if (pool.size() > 0) 
		return pool[chrIDX];

	MinChrom minchrom = minPool[chrIDX];
	Chromosome i;
	i = minchrom;
	return i;
}
	
void ChromPool::ApplyPresentGenotypes(Individual *ind, bool modelLociOnly) {
	if (pool.size() > 0 || minPool.size() > 0) {
		if (ind->GetFather() == NULL) {	
			Chromosome chr1 = GetExpression(ind->ChrID1(chromID));
			Chromosome chr2 = GetExpression(ind->ChrID2(chromID));
			ind->SetChromosomalData(chromID, chr1, chr2);
			//This only does anything if we don't have complete data
			ind->ExpandModelLoci(chromID, modelLoci);
/*chr1.ShowGenotypes(25, ' ');
cout<<"\t";
chr2.ShowGenotypes(25, ' ');
cout<<"\n";
*/

		}
		else 
			ind->ExpandModelLoci(chromID, modelLoci);
	}
}
void ChromPool::PopulatePool(Individual* ind) {
	ind->PopulateChromosomePool(chromID, pool);
}
void ChromPool::DrawIndividual(Individual* ind) {
	ApplyChromosomeIDs(ind);
	ApplyPresentGenotypes(ind, true);
}

void ChromPool::ResolveGenotypes(Individual *ind) {
	if (ind->NeedsCompleteGenotypeInformation())
		ind->SetChromosomalData(chromID, pool[ind->ChrID1(chromID)], pool[ind->ChrID2(chromID)]);
	ind->ResolveGenotypes(chromID, modelLoci);
}

// Clears the chromosomal pool and individuals in this population
// Arg: none
// Ret: none
void ChromPool::Clear(){
	usedIndividuals.clear();
	pool.clear();
}


void ChromPool::BlockDefinition::Randomize() {
	minSnpCount=Utility::Random::globalGenerator.lrand(2, 10);
	maxSnpCount=minSnpCount+Utility::Random::globalGenerator.lrand(2, 45);

	maxBlckMap=(float)Utility::Random::globalGenerator.lrand(1,800)/100000.0;
	minBlckMap=maxBlckMap/(float)Utility::Random::globalGenerator.lrand(10,1000);
	if (minBlckMap < 0.000001)
		minBlckMap = 0.000001;
	
	maxSnpMap=maxBlckMap/(float)Utility::Random::globalGenerator.lrand(1,100);
	minSnpMap=maxSnpMap/(float)Utility::Random::globalGenerator.lrand(10,10000);
	if (minSnpMap < 0.000001)
		minSnpMap = 0.000001;

	frequency=(float)(Utility::Random::globalGenerator.lrand(10,35))/100.0;

	WriteConfiguration(cout);
	cout<<"\n";
}


}
