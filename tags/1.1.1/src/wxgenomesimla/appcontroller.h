//
// C++ Interface: appcontroller
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIAPPCONTROLLER_H
#define GENOMESIM_GUIAPPCONTROLLER_H
#include <string>
#include "appconfig.h"
#include "simulation/growthrate.h"
#include "simulation/ldpngcomponent.h"
#include "simulation/ldplotter.h"
#include "utility/executionlog.h"
#include <pthread.h>
#include <wx/filename.h>
#include "locusmanager.h"
#include "utility/types.h"
namespace GenomeSIM {

namespace GUI {

using namespace std;
using namespace Utility;



/**
 *	@brief Controller object for the application
 * 	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class AppController{
public:
    AppController();
    ~AppController();

	void Refresh();

	static StringArray gaSettings;
	/**
	 * @brief initializes the user's default content- copying files if needed
	 */
	void InitUserResources();

	void ClearLog();
class SimParameters {
public:
	SimParameters();
	~SimParameters();
	
	void DumpToFile();

	void SetConfigFilename(const char *cfg, bool doLoad = true);
	string GetConfigFilename();
	bool WriteConfiguration();
	bool WriteConfiguration(const char *config);
	AppConfig *GetConfiguration() { return configuration; }

	GrowthRate *GetGrowthRate();
	void SetGrowthRate(const char *growthCfg);
	
	uint GetRandomSeed();
	void SetRandomSeed(uint seed);
	uint GetDropCount();
	void SetDropCount(uint count);
	uint GetDropFreq();
	void SetDropFreq(uint freq);
	uint GetDropInit();
	void SetDropInit(uint init);
	uint GetSimultaneousChroms();
	void SetSimultaneousChroms(uint chroms);
	uint GetThreadsPerChrom();
	void SetThreadsPerChrom(uint threads);
	PenetranceModel *GetDiseaseModel();
	void SetDiseaseModel(PenetranceModel* model);
	void DefineDiseaseModel(const char *modelDetails);

	ChromPool *AddChromosomePool(uint blockCount, ChromPool::BlockDefinition &defBlock, const char *label);
	ChromPool *LoadChromosomePool(const char *locSource, const char *label);
	void DeleteChromosomePool(uint idx);
	void ClearLocusSelections();
	void AddLocusSelector(LocusSelection &loc);
	LocusSelectionMap& GetLocusSelections();
	ChromPool::BlockDefinition *GetDefaultBlock();
	void SetDefaultBlock(ChromPool::BlockDefinition& defBlock);


	void SetFontFilename(const char *filename);
	string GetFontFilename();

	void SetMaxSnpDistance(uint distance);
	uint GetMaxSnpDistance();

	void SetStyleSheetFilename(const char *filename);
	string GetStyleSheetFilename();

	void SetReportBufferSize(uint buffSize);
	uint GetReportBufferSize();

	void SetTargetPopulation(uint targetPop);
	uint GetTargetPopulation();

	void WritePairwiseLD(bool doWrite);
	bool WritePairwiseLD();

	void PlotDPrime(bool doPlot);
	bool PlotDPrime();

	void PlotRSquared(bool doPlot);
	bool PlotRSquared();

	void FastLD(bool isFast);
	bool FastLD();

	void SetProjectName(const char *projectName);
	string GetProjectName();

	uint GetStartingGeneration();
	void SetStartingGeneration(uint gen);

	string GetCwd();


	
	
	AppConfig *configuration;
	std::string config;
	std::string cwd;
};

	void CalculateAlleleFrequencies();

	
	/**
	 * @brief Summarize the simulation settings
	 */
	void SummarizeSimulation(ostream &log);

	void SetConfigurationFilename(const char *filename, bool doLoad = true);

	/**
	 * @brief Prepare the simulation to run
	 */
	bool InitializeExecution(uint startGeneration, ostream& log);

	size_t InitializeExecution(ostream& log);
	
	/**
	 * @brief Run simulation for generationCount generations
	 */
	void RunSimulation(uint generationCount, bool doLoadFirst);

	size_t RunSimulation(bool doLoadFirst);
	void PerformAnalysis();
	string OpenSummaryReport(ofstream& os);

	void WriteDatasets();

	SimParameters parameters;
	bool SimIsRunning(); 
	
	float GetPercentCompleted() { return percCompleted; }

	//Return a string containing output from simulation since last read
	string GetSimOutput();

	/**
	 * @brief Returns the current generation of the simulation of pool, pool (threadsafe)
	 */
	uint GetCurGeneration();

	/**
	 * @brief Returns the size of the pool, pool (threadsafe)
	 */
	uint GetCurPopulation();

	/**
	 * @brief Flags the simulation to halt as soon as it can
	 */
	void HaltExecution();

	void GetSessionEntries(vector<ExecutionLog::LogEntry>& entries);
	

	ExecutionLog::RunType *GetProjectEntries();

	ExecutionLog::RunType *GetProjectEntries(const char *project);
//	bool GetProjectEntries(vector<ExecutionLog::LogEntry>& entries);

	bool GetProjectList(vector<string>& entries);

	void InitExecLog();

	void AppendActiveEntry(ExecutionLog::LogEntry &entry);

	void ClearCurrentLogEntries();

	void LockThread();
	void UnlockThead();
	/**
	 * @brief Loads chromosomes (if possible) into the vector
	 */
	bool LoadFileBasedChroms(vector<FileBasedChromosome *> &chroms);
	void SummarizeDiseaseModel(ostream &os, vector<Locus*>&Loci);

	void InitExecution();

	bool InitializationComplete() { return initializationComplete; }

	void UpdateCWD(const char *newCWD);
protected:
	static void *ThreadedInitExecution(void *controller);
	static void *ThreadedSimulation(void *simDetails);
	static void *ThreadedAnalysis(void *controller);
	static void *ThreadedDatasetGeneration(void *controller);
	void GenerateDatasets();


	void AnalyzeCurrentPool();

	std::string config;
	ExecutionLog *execLog;
	PoolManager *pools;	
	bool simIsRunning;
	stringstream simOutput;
	pthread_t simThread;
	pthread_mutex_t simLock;
	bool continueRunning;
	bool initializationComplete;
	vector<ExecutionLog::LogEntry> curEntries;
	float percCompleted;
};

inline
void AppController::SimParameters::SetConfigFilename(const char *cfg, bool doLoad /*=true*/) {
	if (doLoad) {
		if (configuration)
			delete configuration;

		configuration = new AppConfig();
		configuration->LoadSettings(cfg);

		wxFileName filename(cfg);
		if (!wxSetWorkingDirectory(filename.GetPath()))
			cout<<"Issues were encountered trying to set the current working directory to: "<<filename.GetPath()<<"\n";
	}
	else
		configuration->SetAppConfigurationFilename(cfg);
	config = cfg;
}

inline
string AppController::SimParameters::GetConfigFilename() {
	return config;
}

inline
string AppController::SimParameters::GetCwd() {
	return cwd;
}

inline
bool AppController::SimIsRunning() { 	
	return simIsRunning; 
}

}

}

#endif
