//
// C++ Interface: genomesimlawrapper
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIMLAWRAPPER_H
#define GENOMESIMLAWRAPPER_H

#include <pthread.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

pthread_mutex_t paramListMutex = PTHREAD_MUTEX_INITIALIZER;

//I just do it this way, so that pthreads can be compiled away, if you don't want threading....we aren't doing that here, though:O
#define LOCK_PARAMLIST pthread_mutex_lock(&paramListMutex);
#define UNLOCK_PARAMLIST pthread_mutex_unlock(&paramListMutex);



struct params {
	int initPop;
	int carryPop;
	int m;
	float growth;
	float t;
	int seed;
	string offspring;
	string bw;
	int idx;
	int founders;
	float distortion;
	int repeats;
	
	//int max[] = {3, 3, 6, 3, 3};
	
	params(int initPop, int carryPop, int m, float growth, float t, int seed, string bw, int idx, int founders=0, float distortion=0.0, int repeats=0) : initPop(initPop), carryPop(carryPop), m(m), growth(growth), t(t), seed(seed), bw(bw), idx(idx), founders(founders), distortion(distortion), repeats(repeats)
	{
		char cmd[2048];
		
		if (founders>0)
			sprintf(cmd, "gsWrapper param%d %d %d %d %f %f %s %d %d %f %d", this->idx, initPop, carryPop, m, growth, t, bw.c_str(), seed, founders, distortion, repeats);
		else 
			sprintf(cmd, "gsWrapper param%d %d %d %d %f %f %s %d", this->idx, initPop, carryPop, m, growth, t, bw.c_str(), seed);
		cout<<"******* "<<cmd<<"\n";
	}
	
	params() : initPop(0), carryPop(0), m(0), growth(0.0), t(0.0), seed(0), offspring(""), bw(""), idx(0), founders(0), distortion(0.0), repeats(0) { }

	params(const params &other) : initPop(other.initPop), carryPop(other.carryPop), m(other.m), growth(other.growth), t(other.t), seed(other.seed), offspring(other.offspring), bw(other.bw), idx(other.idx), founders(other.founders), distortion(other.distortion), repeats(other.repeats) { }

	void RunGenomeSIMLA() {
		char cmd[2048];
		if (founders>0) {
			sprintf(cmd, "gsWrapper param%d %d %d %d %f %f %s %d %d %f %d", this->idx, initPop, carryPop, m, growth, t, bw.c_str(), seed, founders, distortion, repeats);
			cout<<"<"<<strlen(cmd)<<">\n";			
		}
		else 
			sprintf(cmd, "gsWrapper param%d %d %d %d %f %f %s %d", this->idx, initPop, carryPop, m, growth, t, bw.c_str(), seed);
		LOCK_PARAMLIST;
		cout<<" -- Attempting to call '"<<cmd<<"'\n";
		UNLOCK_PARAMLIST;
		int rv = system(cmd);
		LOCK_PARAMLIST;
		cout<<"-- "<<cmd<<" completed\n";
		UNLOCK_PARAMLIST;
	}
};
/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GenomeSimlaWrapper{
public:
    GenomeSimlaWrapper(uint threadCount);
    ~GenomeSimlaWrapper();

	void Start();

	params *GetNextParam();
	static void *RunThread(void *args);

	void LoadParameters(const char *params);

protected:
	vector<params> parameterSettings;
	static uint curIdx;
	static bool UseAdamEve;
	uint threadCount;

	vector<int> initPop;
	vector<int> carryCap;
	vector<int> M;
	vector<float> growth;
	vector<float> t;
	vector<int> seeds;
	vector<string> offspring;
	vector<string> bw;
	vector<int> founderCount;
	vector<float> distortion;
	vector<int> repeats;


};

#endif
