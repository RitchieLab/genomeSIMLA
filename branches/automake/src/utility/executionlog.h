//
// C++ Interface: executionlog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIMEXECUTIONLOG_H
#define GENOMESIMEXECUTIONLOG_H

#include "types.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "exception.h"

namespace Utility {

using namespace std;
	

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ExecutionLog{
public:

//I'm being bad and adding in some genomeSIMLA specific stuff. Those pieces really need to be put into a
//class that inherits from this one and the structures that hold these items should be templated accordingly
struct LogEntry {
	size_t startingGeneration;
	size_t currentGeneration;
	uint randomSeed;
	time_t startTime;
	time_t endTime;
	time_t duration;
	string reportFilename;
	string projectName;
	uint expressionCount;

	//GnomeSIMLA Stuff
	bool sampledReport;
	uint blockCount;
	

	map<string, string> reports;

	LogEntry() : startingGeneration(0), currentGeneration(0), randomSeed(0), reportFilename(""), projectName(""), sampledReport(false), blockCount(0)	{}

	LogEntry(size_t startGen, size_t curGen, uint seed, const char *filename, const char *projectName) : 
			startingGeneration(startGen), currentGeneration(curGen), randomSeed(seed), reportFilename(filename), 
			projectName("") {
		startTime = GetTime();
	}
	
	LogEntry& operator=(const LogEntry& other) {
		startingGeneration 	= other.startingGeneration;
		currentGeneration	= other.currentGeneration; 
		randomSeed			= other.randomSeed;
		startTime			= other.startTime;
		endTime				= other.endTime;
		duration			= other.duration;
		reportFilename		= other.reportFilename;
		projectName			= other.projectName;
		expressionCount		= other.expressionCount;
		sampledReport		= other.sampledReport;
		blockCount			= other.blockCount;
		reports				= other.reports;
		return *this;
	}

	LogEntry(const LogEntry& other) : 
			startingGeneration(other.startingGeneration), currentGeneration(other.currentGeneration), randomSeed(other.randomSeed), 
			startTime(other.startTime), endTime(other.endTime), duration(other.duration), reportFilename(other.reportFilename),
			projectName(other.projectName), expressionCount(other.expressionCount), sampledReport(other.sampledReport), blockCount(other.blockCount),
			reports(other.reports) { }

	void Start(size_t startGen, size_t currGen, uint seed, const char *filename, const char *projectName) {
		startTime = GetTime();
		startingGeneration = startGen;
		randomSeed = seed;
		reportFilename = filename;
		currentGeneration = currGen;
		this->projectName = projectName;
		endTime = duration = 0;
	}

	void Stop(uint expCount) {
		endTime = GetTime();
		duration = endTime - startTime;
		expressionCount = expCount;
	}

	bool IsValid() {
		return startingGeneration != 0 || currentGeneration != 0;
	}

	time_t GetTime() {
		return time(NULL);
	}

	friend std::istream &operator>>(std::istream &is, LogEntry &entry);
	friend std::ostream &operator<<(std::ostream &os, LogEntry &entry);


};

	/**
	 * @brief Map containing each of the runs for a given project
	 * The key associated with the map is the generation. It is 
	 * assumed that the most recent entry will be retained. 
	 */
	typedef map<int, LogEntry> RunType;

	/**
	 * @brief Map containing the various projects associate with a 
	 * single configuration. The key is the project name. 
	 */
	typedef map<string, RunType> EntryType;
	typedef pair<string, LogEntry> PairType;

	/**
	 * @brief Maintains the log of executions. 
	 */
    ExecutionLog(const char *configFilename);

    ~ExecutionLog();

	void Initialize();

	/**
	 * @brief Initialize this set of runs
	 * @note If startingGeneration > 0, we will retain all of the previous entries. Otherwise, they will be removed
	 */
	void InitProject(uint seed, uint startingGeneration, const char *projectName);

	/**
	 * @brief Initiates a new run. 	
	 * @param generation The generation of the last generation to be simulated
	 * @note We will record the start time, generation, starting generation (which might have been 0), random seed
	 */
	void StartRun(uint generation);

	/**
	 * @brief Indicates that the run is completed
	 * @param expressionCount The size of each chromosome pool when the run is compeleted
	 * @note We will record the starting time and the expressionCount
	 */
	LogEntry EndRun(uint expressionCount, const char *reportFilename, const char *index);


	LogEntry EndRun(uint generation, uint expCount, const char *reportFilename, const char *index);

	/**
	 * @brief Save the log out to file
	 */
	void SaveFile();

	/**
	 * @brief Remove all entries associated with project, project.
	 */
	void ClearEntries(const char *project);

	void ClearEntries();

	/**
	 * @brief Append a report of type, reportType, to the current log entry
	 */
	void AppendReport(const char *reportType, const char *filename);

	/**
	 * @brief Returns a vector containing each entry associated with the project, project.
	 * @return T/F depending on whether or not anything was found for that project name
	 */
	RunType *GetProjectEntries(const char *project);

	/**
	 * @brief Returns the list of projects found in the log
	 */
	bool GetProjectList(vector<string> &projects);

	void UpdateEntry(const char *project, ExecutionLog::LogEntry &entry);

protected:
	void LoadFile();

	string filename;					///<This is the actual log filename (based on the configuration filename)
	string projectName;					///<This is the name of the project name to be used. We cache it to pass it on to new entries
	uint seed;							///<The random seed being used
	uint startingGeneration;			///<The first generation of this run(s)
	LogEntry currentEntry;				///<The current entry we are working on 
	ofstream log;						///<the file the log is being written to
	
	/**
	 * Each project's runs associated with the project name
	 */
	EntryType entries;	
	//vector<LogEntry> entries;			///<The various entries we are using
};

	
		
	
	
}


#endif
