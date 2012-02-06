//
// C++ Implementation: executionlog
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "executionlog.h"
#include "strings.h"


namespace Utility {

ExecutionLog::ExecutionLog(const char *configFilename) : seed(0), startingGeneration(0) {
	stringstream ss;
	ss<<configFilename<<".log";
	this->filename = ss.str();
}


ExecutionLog::~ExecutionLog()	{
	log.close();
}




/**
 * This is only required if something has been deleted
 */
void ExecutionLog::SaveFile() {
	//Reset to the beginning of the file stream
	log.close();
	log.open(this->filename.c_str(), ios::out|ios::trunc);

	EntryType::iterator projectItr = entries.begin();
	EntryType::iterator projectEnd = entries.end();

	while (projectItr != projectEnd) {
		RunType::iterator itr = projectItr->second.begin();
		RunType::iterator end = projectItr->second.end();
		while (itr != end) {
			log<<itr->second<<"\n";
			itr++;
			log.flush();
		}
		projectItr++;
	}
}
	
void ExecutionLog::LoadFile() {
	//First, we want to open the file (if it exists) and read in all entries
	ifstream prevLog(filename.c_str());

//	cout<<"Opening File: "<<filename<<"\n";

	//If it's there and has data, let's save it
	if (prevLog.good()) {
		while (!prevLog.eof()) {
			LogEntry newEntry;
			char line[4096];
			prevLog.getline(line, 4096);
			
			if (strlen(line) > 0){
//				cout<<"strlen(line): "<<strlen(line)<<"\n";
				stringstream ss(line);
				ss>>newEntry;
//				cout<<" > "<<newEntry.projectName<<" : "<<newEntry<<"\n";
				if (newEntry.IsValid())
				(entries[newEntry.projectName])[newEntry.currentGeneration] = newEntry;
			
				if (newEntry.projectName == "") {
//					cout<<"Throwing FileIOError exception\n";
					throw Exception::BadFileVersion(filename.c_str());
				}
			}
		}
	}
}


void ExecutionLog::ClearEntries(const char *projectName) {
	entries.erase(projectName);
}

void ExecutionLog::ClearEntries() {
	entries.clear();
}


void ExecutionLog::Initialize() {
	LoadFile();
	log.open(this->filename.c_str(), ios::out|ios::app);
}

void ExecutionLog::InitProject(uint seed, uint startingGeneration, const char *projectName) {
	this->seed = seed;
	this->startingGeneration = startingGeneration;
	this->projectName = projectName; 	//ExtractFilename(projectName);
	
	if (startingGeneration == 0) {
		//We want to clear any previous entries
		log.close();
		log.open(filename.c_str(), ios::out);
		ClearEntries(this->projectName.c_str());
		SaveFile();
	}

//	cout<<"Exec::Log : Init("<<startingGeneration<<")\n";

}

ExecutionLog::RunType *ExecutionLog::GetProjectEntries(const char *project) {
	RunType *prjEntries = NULL;

	EntryType::iterator itr = entries.find(project);
	if (itr != entries.end()) 	
		prjEntries = &itr->second;

	return prjEntries;
}

bool ExecutionLog::GetProjectList(vector<string> &projects) {
	EntryType::iterator itr = entries.begin();
	EntryType::iterator end = entries.end();

	string lastKey = "";
	while (itr != end) {
//		cout<<"-- Project "<<itr->first<<" : "<<itr->second<<"\n";
//		if (itr->first != lastKey) {
			lastKey = itr->first;
			projects.push_back(lastKey);
//		}
		itr++;
	}
	return projects.size() > 0;
}


void ExecutionLog::StartRun(uint generation) {
//	cout<<"Exec::Log : Start Run ("<<generation<<")]\n";
	currentEntry.Start(startingGeneration, generation, seed, filename.c_str(), projectName.c_str());
}

ExecutionLog::LogEntry ExecutionLog::EndRun(uint generation, uint expCount, const char *reportFilename, const char *index) {
	startingGeneration = generation;
	
	return EndRun(expCount, reportFilename, index);
}

ExecutionLog::LogEntry ExecutionLog::EndRun(uint expressionCount, const char *reportFilename, const char *index) {
	currentEntry.Stop(expressionCount);
	currentEntry.reportFilename = reportFilename;
	currentEntry.reports[index] = reportFilename;
	(entries[projectName])[currentEntry.currentGeneration]=currentEntry;
	log<<currentEntry<<"\n";

	log.flush();
//	cout<<"Exec::Log : End Run ("<<currentEntry.currentGeneration<<")\n";
	return currentEntry;
}

void ExecutionLog::AppendReport(const char *reportType, const char *reportFilename) {
	currentEntry.reports[reportType] = reportFilename;
}


void ExecutionLog::UpdateEntry(const char *project, ExecutionLog::LogEntry& entry){
	RunType &projectEntries = entries[project];
	projectEntries[entry.currentGeneration] = entry;
}


std::ostream &operator<<(std::ostream &os, ExecutionLog::LogEntry &entry) {
	os<<entry.startingGeneration<<" "
		<<entry.currentGeneration<<" "
		<<entry.randomSeed<<" "
		<<entry.startTime<<" "
		<<entry.endTime<<" "
		<<entry.duration<<" "
		<<entry.expressionCount<<" \""
		<<entry.projectName<<"\" "
		<<entry.sampledReport<<" "
		<<entry.blockCount<<" ";
	map<string,string>::iterator itr = entry.reports.begin();
	map<string,string>::iterator end = entry.reports.end();
		
	while (itr != end) {
		os<<itr->first<<" \""<<itr->second<<"\" ";
		itr++;
	}
	return os;
}

std::istream &operator>>(std::istream &is, ExecutionLog::LogEntry &entry) {
	is>>entry.startingGeneration
		>>entry.currentGeneration
		>>entry.randomSeed
		>>entry.startTime
		>>entry.endTime
		>>entry.duration
		>>entry.expressionCount;
	entry.projectName = ParseFilename(is);
	is>>entry.sampledReport		
		>>entry.blockCount;
	while (!is.eof()){
		string key="", value="";
		is>>key;
		value = ParseFilename(is);

		if (key.length() > 0 && value.length() > 0)
			entry.reports[key] = value;
		else
			break;
	}
	entry.reportFilename = entry.reports["Index"];
	return is;
}

}

