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
	distribution[snp->GetLabelCount() - 1].push_back(DistContents(snp->GetLabel(), value, id));
	MPIUNLOCK;
}

uint NTest::GetTestCount(uint modelSize) {
	assert(modelSize < maxModelSize);
	return distribution[modelSize].size();
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
				*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3)<<setw(columnWidth)<<i;
				*os<<setw(columnWidth)<<distribution[h][i].testID;
				*os<<setw(columnWidth)<<distribution[h][i].modelLabel<<setw(columnWidth+6)<<distribution[h][i].value<<"\n";
			}
		}
		*os<<"\n";
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

	while (distribution[idx][i++].value > value && i<count);
	return (float)i/(float)count;
}


void NTest::Sort() {
	for (uint i=0; i<maxModelSize; i++) 
		sort(distribution[i].begin(), distribution[i].end());

	isSorted=true;
}




/************************************************ Omnibus Distribution *********************/
Omnibus::Omnibus(int seed) : PTestDistribution("", seed) { }
Omnibus::~Omnibus() { }
void Omnibus::AppendTest(float value, uint modelSize, const char *modelLabel, int testID) {
	isSorted=false;
	int id=testID;
	if (testID < 0)
		id = distribution.size();
	distribution.push_back(DistContents(modelLabel, value, testID));
}

void Omnibus::AppendTest(float value, SnpAligned *snp, int testID) {
	int id=testID;
	if (testID < 0)
		id = distribution.size();
	isSorted=false;
	distribution.push_back(DistContents(snp->GetLabel(), value, testID));
}

float Omnibus::GetPValue(float value, uint modelSize) {
	if (!isSorted) 
		Sort();

	uint i=0;
	uint count=distribution.size();

	while (distribution[i++].value > value && i<count);
	return (float)i/(float)count;
}

void Omnibus::Sort() {
	sort(distribution.begin(), distribution.end());
	isSorted = true;
}
uint Omnibus::GetTestCount(uint modelSize) {
	return distribution.size();
}

}

}
