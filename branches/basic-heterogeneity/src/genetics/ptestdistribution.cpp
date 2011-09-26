//
// C++ Implementation: ptestdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ptestdistribution.h"
#include <iomanip>



#ifdef USE_MPI
#include <pthread.h>
pthread_mutex_t mpi_dist_lock  = PTHREAD_MUTEX_INITIALIZER;
#define MPILOCK pthread_mutex_lock(&mpi_dist_lock)
#define MPIUNLOCK pthread_mutex_unlock(&mpi_dist_lock)

#else
//We don't want these to mean anything if we arne't using MPI
#define MPILOCK
#define MPIUNLOCK

#endif


namespace Genetics {

namespace Distribution {


/******************************   N Test Distribution *****************************************/
NTest::NTest(uint modelSize, int seed, char *label) : PTestDistribution(label, seed), maxModelSize(modelSize) {
	distribution = new vector<DistContents>[modelSize];
}

NTest::~NTest() {
	if (distribution)
		delete[] distribution;
}

void NTest::AppendTest(float value, uint modelSize, const char *modelLbl, int testID) {
	MPILOCK;
	int id=testID;
	if (testID < 0)
		id = distribution[modelSize].size();

	if (maxModelSize <= modelSize)
		cout<<"Ugh! We can't fit a model of size: "<<modelSize<<" into the "<<label<<" distribution! It's only for models of size: "<<maxModelSize<<" or smaller!\n";
	assert(modelSize < maxModelSize);
	//cout<<"Appending Test"<<testID<<": "<<modelLbl<<" - "<<value<<" - "<<id<<"\n";

	distribution[modelSize].push_back(DistContents(modelLbl, value, id));
	MPIUNLOCK;
}

void NTest::AppendTest(float value, SnpAligned *snp, int testID) {
	MPILOCK;
	int id=testID;
	uint modelSize = snp->GetLabelCount();
	if (testID < 0)
		id = distribution[modelSize].size();
	assert(snp->GetLabelCount() <= maxModelSize);
//	cout<<"Appending Test"<<testID<<": "<<snp->GetID()<<" - "<<value<<"\n";
	distribution[snp->GetLabelCount() - 1].push_back(DistContents(snp->GetLabel(), value, id));
	MPIUNLOCK;
}

uint NTest::GetTestCount(uint modelSize) {
	assert(modelSize < maxModelSize);
	return distribution[modelSize].size();
}

void Omnibus::LoadDistribution(istream *os) {
	//Let's dump anything that might be left behind
	distribution.clear();

	char line[4096];
	
	//We can ignore the first 3 lines, since they are associated with header information
	for (uint i = 0; i<3; i++) {
		os->getline(line, 4096);
	}

	int idx;
	//For this distribution, we know how many there will be and they have a predefined order.
	while (!os->eof()) {
		DistContents c;

		string val;
		
		*os>>idx;
		if (!os->eof() && idx > 0){

			*os>>c.testID>>c.modelLabel>>val;
			c.value = 0.0;
			if (strcmp(val.c_str(), "inf") == 0) {
				c.value=1000.0/c.value;
			}
			else if (strcmp(val.c_str(), "nan") == 0) {
				c.value = nanf("char-sequence");
			}
			else {
				stringstream ss(val);
				//cout<<"Trying to extract "<<val<<" to the structure\n";
				ss>>c.value;
			}
			
			distribution.push_back(c);
//			distribution[idx]=c;
		}
	}

	
}
void Omnibus::DumpDistribution(ostream* os) {
	Sort();
	if (os) {
		uint columnWidth=8;
		*os<<label<<" distribution\n";
		uint count=distribution.size();

		*os<<setw(columnWidth)<<" "<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Model"<<setw(columnWidth+6)<<label<<endl;
		*os<<setw(columnWidth)<<"ID"<<setw(columnWidth)<<"Number"<<setw(columnWidth)<<" ID "<<setw(columnWidth+6)<<"P Value"<<endl;
		for (uint i=0; i<count; i++) {
			*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3)<<setw(columnWidth)<<i + 1;
			*os<<setw(columnWidth)<<distribution[i].testID + 1;
			*os<<setw(columnWidth)<<distribution[i].modelLabel<<setw(columnWidth+6)<<distribution[i].value<<"\n";
		}
		*os<<"\n";
	}
}

void NTest::LoadDistribution(istream *os) {
	//Let's dump anything that might be left behind
	for (uint i=0; i<maxModelSize; i++) 
		distribution[i].clear();

	char line[4096];
	//First line is garbage
	os->getline(line, 4096);

	for (uint i=0; i<maxModelSize; i++)  {
		//Next three lines are header details
		os->getline(line, 4096);
		os->getline(line, 4096);
		os->getline(line, 4096);

		int idx = 1;
		//For this distribution, we know how many there will be and they have a predefined order.
		while (!os->eof() && idx > 0) {
			DistContents c;
			idx=0;
			*os>>idx;
			if (!os->eof() && idx > 0){
				*os>>c.testID>>c.modelLabel>>c.value;
	
				distribution[i].push_back(c);
			}
		}
	}
}


void NTest::DumpDistribution(ostream* os) {
	Sort();
	if (os) {
		uint columnWidth=8;
		*os<<label<<" distribution\n";
		for (uint h=0; h<maxModelSize; h++) {
			uint count=distribution[h].size();
			*os<<h+1<<" SNP models\n";
			*os<<setw(columnWidth)<<" "<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Model"<<setw(columnWidth+6)<<label<<endl;
			*os<<setw(columnWidth)<<"ID"<<setw(columnWidth)<<"Number"<<setw(columnWidth)<<" ID "<<setw(columnWidth+6)<<"P Value"<<endl;
			for (uint i=0; i<count; i++) {
				*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3)<<setw(columnWidth)<<i + 1;
				*os<<setw(columnWidth)<<distribution[h][i].testID + 1;
				*os<<setw(columnWidth)<<distribution[h][i].modelLabel<<setw(columnWidth+6)<<distribution[h][i].value<<"\n";
			}
		}
		*os<<" \n";
	}
}


float NTest::GetPValue(float value, uint modelSize) {
	if (modelSize > maxModelSize) {
		cout<<"Internal error. Distribution call made with an inappropriate size of "<<modelSize<<"/"<<maxModelSize<<"\n";
		abort();
	}
	if (!isSorted) 
		Sort();

	uint i=0;
	uint idx=modelSize;
	uint count=distribution[idx].size();

	while (i<count && distribution[idx][i++].value > value);
	return (float)i/(float)count;
}


void NTest::Sort() {
	for (uint i=0; i<maxModelSize; i++) 
		sort(distribution[i].begin(), distribution[i].end());

	isSorted=true;
}




/************************************************ Omnibus Distribution *********************/
Omnibus::Omnibus(int testCount, int seed, const char *label, float defInitValue /*= 0.0*/) : PTestDistribution(label, seed), distribution(testCount, DistContents("", defInitValue, 0)), lastID(0), testCount(testCount) {

}

Omnibus::~Omnibus() { }
void Omnibus::AppendTest(float value, uint modelSize, const char *modelLabel, int testID) {
	uint id=testID;
	if (testID < 0)
		id = lastID++;

	assert((uint)testID<distribution.size());
		
	isSorted=false;

	if (distribution[testID].value < value)
		distribution[testID] = DistContents(modelLabel, value, id);
}

void Omnibus::AppendTest(float value, SnpAligned *snp, int testID) {
	uint id=testID;
	if (testID < 0)
		id = lastID++;
	isSorted=false;
	
	assert((uint)testID<distribution.size());

	if (distribution[testID].value < value)
		distribution[testID] = DistContents(snp->GetLabel(), value, id);
}



float Omnibus::GetPValue(float value, uint modelSize) {
	if (!isSorted) 
		Sort();

	uint i=0;
	uint count=distribution.size();

	if (count == 0)
		return 0.0;
	while (i<count && distribution[i++].value > value );
	return (float)i/(float)count;
}

void Omnibus::Sort() {
	sort(distribution.begin(), distribution.end());
	isSorted = true;
}
uint Omnibus::GetTestCount(uint modelSize) {
	return distribution.size();
}


void InvertedOmnibus::AppendTest(float value, uint modelSize, const char *modelLabel, int testID) {
	uint id=testID;
		if (testID < 0)
	id = lastID++;

	assert((uint)testID<distribution.size());
		
	isSorted=false;

	if (value < distribution[testID].value )
		distribution[testID] = DistContents(modelLabel, value, id);
}
void InvertedOmnibus::AppendTest(float value, SnpAligned *snp, int testID) {
	uint id=testID;
	if (testID < 0)
		id = lastID++;
	isSorted=false;
	
	assert((uint)testID>distribution.size());

	if (value < distribution[testID].value)
		distribution[testID] = DistContents(snp->GetLabel(), value, id);
}


class gt {
public:
	bool operator()(const DistContents& l, const DistContents& r) { return l.value<r.value; }
};

void InvertedOmnibus::Sort() {
	sort(distribution.begin(), distribution.end(), gt());
	isSorted = true;
}
float InvertedOmnibus::GetPValue(float value, uint modelSize) {
	if (!isSorted) 
		Sort();

	uint i=0;
	uint count=distribution.size();

	if (count == 0)
		return 0.0;
	while (i<count && distribution[i++].value < value );
	return (float)i/(float)count;
}

}

}
