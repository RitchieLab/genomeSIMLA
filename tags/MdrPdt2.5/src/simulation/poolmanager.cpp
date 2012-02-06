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

namespace Simulation {



string PoolManager::pathToJava 			= "/usr/local/bin/java";
string PoolManager::pathToHaploview 	= "/usr/local/bin/Haploview.jar";
string PoolManager::haploviewSettings 	= "-png -n";
string PoolManager::haploviewOverviewSettings = "-compressedpng -n";
string PoolManager::javaSettings		= "-Xms1000m -Xmx2000m";
int PoolManager::haploWindowSize		= 199;
int PoolManager::haploWindowStride		= 130;
bool PoolManager::generateOverviewLD	= true;

//bool PoolManager::useDistanceMaps			= true;
//bool PoolManager::distanceMappingFunction	= NULL;


using namespace PopulationGrowth;

PoolManager::PoolManager() : currGeneration(0), poolIsPopulated(false)
{
}


PoolManager::~PoolManager() {
	uint count = pools.size();
	for (uint i=0; i<count; i++) 
		delete pools[i];
}


void PoolManager::AddPool(ChromPool *pool) {
	pools.push_back(pool);
}

string PoolManager::GenerateMarkerInfoFilename(uint pos, uint segmentID, const char *prj) {
	stringstream ss;
	ss<<prj<<"."<<pos<<"."<<segmentID<<".info";
	return ss.str();
}

string PoolManager::GenerateMarkerInfoFilename(uint pos, const char *prj) {
	stringstream ss;
	ss<<prj<<"."<<pos<<".info";
	return ss.str();
}

string PoolManager::GenerateLocusFilename(uint gen, uint pos, const char *prj) {
	stringstream ss;
	ss<<prj<<"."<<gen<<"."<<pos<<".loc";
	return ss.str();
}

string PoolManager::GeneratePhasedFilename(uint gen, uint pos, const char *prj) {
	stringstream ss;
	ss<<prj<<"."<<gen<<"."<<pos<<".phased";
	return ss.str();
}

string PoolManager::GeneratePhasedFilename(uint gen, uint pos, const char *prj, int start) {
	stringstream ss;
	ss<<prj<<"."<<gen<<"."<<pos<<"."<<start<<".phased";
	return ss.str();
}

bool PoolManager::VerifyHaploview() {
	cout<<"Haploview verification is under construction\n";
	return true;
}


bool PoolManager::GenerateLDMap(const char *settings, const char *filename, const char *markerinfo, bool leavePhased) {
	//Call haploview on the file we just generated
	bool success = false;
	stringstream cmd;
	if (VerifyHaploview()) {
		cmd<<pathToJava<<" "<<javaSettings<<" -jar "<<pathToHaploview<<" -haps "<<filename<<" "<<settings;
		if (markerinfo)
			cmd<<" -info "<<markerinfo;
		cmd<<"\n";

		cout<<"Executing: "<<cmd.str()<<"\n";

		int rv = system(cmd.str().c_str());
		if (rv)
			cout<<"Unable to execute haploview's LD conversion on "<<filename<<".\n";
		else
			success = true;
	}
	if (!leavePhased) {
		unlink(filename);
	}

	return success;
}

bool PoolManager::GenerateLDMapOverview(const char *filename, const char *markerinfo, bool leavePhased) {
	GenerateLDMap(haploviewOverviewSettings.c_str(), filename, markerinfo, leavePhased);
	return true;
}

bool PoolManager::GenerateLDMap(const char *filename, const char *markerinfo, bool leavePhased) {
	GenerateLDMap(haploviewSettings.c_str(), filename, markerinfo, leavePhased);
	return true;	
}

bool PoolManager::GenerateLDMaps(const char *prj, bool leavePhased) {
	bool success = true;
	uint count = pools.size();
	cout<<"Generating LD Maps for "<<count<<" chromosomes. ";
	if (leavePhased)
		cout<<"Phased data will be left behind. \n";
	else
		cout<<"Phased data will be removed. \n";


	//Work through the various chromosome pools:
	// 1 Produce the overview image, if requested
	// 2 Work through the windows (increment by stride and show windows size number of snps 
	//   per window
	// 3 Produce the locus reports
	// * The phased output is always removed after we haploview generates the png
	// * The main overview pool is left behind if leavePhased is true
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		string filename = GeneratePhasedFilename(currGeneration, i + 1, prj);
	
		ofstream file(filename.c_str(), ios_base::out);
		ch->SaveAsPhased( file );
		
		if (!file.is_open()){
			cout<<"Unable to write to chromosome data file: "<<filename<<". Verify that there is enough disk space and permissions are set properly. No LD maps were generated for this \n";
			success=false;
			break;
		}

		file.close();

		if (generateOverviewLD) {
			string sLocFilename = GenerateMarkerInfoFilename(i, prj);
//			locFilename = sLocFilename.c_str();
//			file.open(locFilename, ios_base::out);
//			ch->WriteMarkerInfo(file, 0, ch->GetLociCount());
			//WriteLocusFile(file, (uint)curr, (uint)(curr+haploWindowSize), ch, *distanceMappingFunction);
//			file.close();
			
			GenerateLDMap(haploviewOverviewSettings.c_str(), filename.c_str(), sLocFilename.c_str(), leavePhased);
		}
		else if (!leavePhased)
			unlink(filename.c_str());
		
	
		int windowStart = 0;				///<The snp idx where the image starts
		int snpCount = ch->GetLociCount();	///<Number of total snps
		int lastSnpWritten = 1;				///<Which was the last snp written
		int segID = 0;						///<The index into the "pieces" of the whole picture
		int windowStartPosition = 0;		///<This is the bp position where the window starts
		DistMappingFn *fn = ChromPool::mappingFn;
		/**
	 	 * Where do we anticipate the next window starting, if there is to be one
	 	 */
		int nextWindowStartPos	= haploWindowStride;

		for (int curr=windowStart; lastSnpWritten<snpCount && lastSnpWritten>0; ) {
			filename = GeneratePhasedFilename(currGeneration, i+1, prj, curr);
			file.open(filename.c_str(), ios_base::out);
			lastSnpWritten = ch->SaveAsPhased(file, curr, haploWindowSize);
			file.close();

			if (lastSnpWritten > 0) {
				const char* locFilename = NULL;
				string sLocFilename;
				if (ch->UseMapInfoFile()) {
					sLocFilename = GenerateMarkerInfoFilename(i, segID++, prj);
					locFilename = sLocFilename.c_str();
					file.open(locFilename, ios_base::out);
					if (fn)
						fn->Reset(windowStartPosition);

					ch->WriteMarkerInfo(file, curr, (curr+haploWindowSize));
					//WriteLocusFile(file, (uint)curr, (uint)(curr+haploWindowSize), ch, *distanceMappingFunction);
	
					file.close();
				}
				GenerateLDMap(haploviewSettings.c_str(), filename.c_str(), locFilename, false);
				curr+=haploWindowStride;
				nextWindowStartPos+=haploWindowStride;

			}
		}

		ch->CalculateAlleleFrequencies();
		filename = GenerateLocusFilename( currGeneration, i, prj);
		file.open(filename.c_str(), ios_base::out);
		if (file.is_open()) {
			ch->SaveLoci( file );
			file.close();
		}
		else {
			cout<<"Unable to write to locus log: "<<filename<<". Verify that there is enough disk space and permissions are set properly.\n";
			success=false;
			break;
		}
		
	}
	return success;
}


bool PoolManager::Dump(const char *prj) {
	bool success = true;
	uint count = pools.size();
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		string filename = GeneratePhasedFilename(currGeneration, i + 1, prj);
		ofstream file(filename.c_str(), ios_base::out);
		if (file.is_open()) {
			ch->SaveAsPhased( file );
			file.close();
		}
		else {
			cout<<"Unable to write to chromosome data file: "<<filename<<". Verify that there is enough disk space and permissions are set properly.\n";
			success=false;
			break;
		}

		ch->CalculateAlleleFrequencies();
		filename = GenerateLocusFilename( currGeneration, i+1, prj);
		file.open(filename.c_str(), ios_base::out);
		if (file.is_open()) {
			ch->SaveLoci( file );
			file.close();
		}
		else {
			cout<<"Unable to write to locus log: "<<filename<<". Verify that there is enough disk space and permissions are set properly.\n";
			success=false;
			break;
		}
		
	}
	return success;
}			

ChromPool *PoolManager::AddChromosome(const char *locFilename) {
	ChromPool *ch = new ChromPool(pools.size(), locFilename);
	pools.push_back(ch);
	return ch;
}

ChromPool *PoolManager::AddChromosome(uint blockCount, float minR, float maxR, ChromPool::BlockDefinition &defaultBlock) {
	ChromPool *ch = new ChromPool(pools.size(), blockCount, minR, maxR);
	ch->DefineDefaultBlock( defaultBlock );
	pools.push_back(ch);
	return ch;
}

	

bool PoolManager::Load( const char *project, uint generation ) {
	bool success = true;

	for (uint i=0; i<pools.size(); i++) {
		ChromPool *ch = pools[i];
		string filename = GenerateLocusFilename(generation, i + 1, project);
		ifstream file(filename.c_str(), ios_base::in );
		if (file.is_open()) {
			ch->LoadLoci( file );
			file.close();
			cout<<"Locus log loaded: \n";
			//ch->SaveLoci( cout );
		}
		else {
			cout<<"Unable to open loci data file: "<<filename<<". Unable to continue\n";
			success=false;	
			break;
		}


		filename = GeneratePhasedFilename(generation, i + 1, project);
		file.open(filename.c_str(), ios_base::in);
		if (file.is_open()) {
			ch->LoadPhased( file, generation );
			file.close();
		}
		else {
			cout<<"Unable to open chromosome data file: "<<filename<<". Verify that the generation data exists to be loaded\n";
			success=false;
			break;
		}
	}
	
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


bool PoolManager::GetAlleleFrequency(uint chrID, uint locId, float &al1, float &al2) {
	bool success = false;
	
	if (chrID < pools.size() ) 
		success = pools[chrID]->GetAlleleFrequencies(locId, al1, al2);
	return success;
}

bool PoolManager::AdvanceGenerations(uint moreGenerations, GrowthRate *f) {
	bool success = true;
	uint count = pools.size();
	for (uint i=0; i<count; i++) {
		ChromPool *ch = pools[i];
		ch->AdvanceGenerations( moreGenerations, f );
	}
	currGeneration+=moreGenerations;
	return success;
}


void PoolManager::CreateInitialPopulation(uint populationCount) {
	uint count = pools.size();
	for (uint i=0; i<count; i++) {
		ChromPool *pool = pools[i];
		pool->BuildPool(populationCount);
	}
	poolIsPopulated = true;
}


void PoolManager::GenerateReport(ostream &os, uint headerWidth) {
	uint count = pools.size();
	cout<<setw(headerWidth)<<"Chromosome Pool Count: "<<count<<endl;

	for (uint i=0; i<count; i++) {
		cout<<setw(headerWidth)<<"Chromosome: "<<pools[i]->GetID()<<"\n";
		cout<<setw(headerWidth + 5)<<"Loci Count: "<<pools[i]->GetLociCount()<<"\n";
		cout<<setw(headerWidth + 5)<<"Expression Count: "<<pools[i]->GetExpressionCount()<<"\n";
	}
}


bool PoolManager::InitializePools(uint startGen, uint poolSize, const char *prj) {
	uint count = pools.size();

	//Create our own loci if we are starting at the beginning
	if (startGen == 0) {
		for (uint i=0; i<count; i++) {
			ChromPool *pool = pools[i];
			pool->InitLoci();

			cout<<"Pool: "<<pool->GetID()<<" Initialized with "<<pool->GetLociCount()<<"\n";
			string locFilename = GenerateMarkerInfoFilename(i, prj);
			ofstream file(locFilename.c_str(), ios_base::out);
			
			if (pool->UseMapInfoFile())
				pool->WriteMarkerInfo(file);
			//This time, we can afford for the fn to be anonymous, since we won't reuse it
			//WriteLocusFile(file, (uint)0, (uint)pool->GetLociCount(), pool, Kosambi());
			file.close();

		}
		for (uint i=0; i<forcedFrequencies.size(); i++) {
			ForcedAF f=forcedFrequencies[i];
			if (f.chrID >= pools.size() || !pools[f.chrID]->ForceAlleleFrequency(f.locus, f.af1, f.af2))
				cout<<"An error was encountered trying to force allele frequencies: "<<f.chrID<<" "<<f.locus<<" "<<f.af1<<" "<<f.af2<<"\n";
		}
		CreateInitialPopulation( poolSize );
		Dump(prj);
	
		poolIsPopulated = true;

	} 
	else {
		poolIsPopulated = Load(prj, startGen);
	}
	
	return poolIsPopulated;

}


}
