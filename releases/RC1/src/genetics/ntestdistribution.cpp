//
// C++ Implementation: ntestdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ntestdistribution.h"

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

namespace Reporting {

void NTestDistribution::Report(ostream* os) {
	uint columnWidth=8;
	uint arrayCount=scores.size();
	uint count;
	for (uint h=0; h<arrayCount; h++) {
		count=scores[h].size(); 
		*os<<h+1<<" SNP models\n";
		*os<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Model"<<setw(columnWidth)<<""<<endl;
		*os<<setw(columnWidth)<<"Rank"<<setw(columnWidth)<<"Number"<<setw(columnWidth)<<" ID "<<setw(columnWidth)<<"Score"<<endl;
		for (uint i=0; i<count; i++) { 
			*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3)<<setw(columnWidth)<<scores[h][i].id;
			*os<<setw(columnWidth)<<scores[h][i].label<<setw(columnWidth)<<scores[h][i].score<<"\n";
		}
	}

}

void NTestDistribution::Report() {
	ostream *os = NULL;
	if (log && log->GetStream())
		os = log->GetStream();
	if (os) {
		uint columnWidth=8;
		uint arrayCount=scores.size();
		uint count;
		for (uint h=0; h<arrayCount; h++) {
			count=scores[h].size(); 

			*os<<h+1<<" SNP models\n";
			*os<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Model"<<setw(columnWidth)<<""<<endl;
			*os<<setw(columnWidth)<<"Rank"<<setw(columnWidth)<<"Number"<<setw(columnWidth)<<" ID "<<setw(columnWidth)<<"Score"<<endl;

			for (uint i=0; i<count; i++) { 
				*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3)<<setw(columnWidth)<<i;
				*os<<setw(columnWidth)<<scores[h][i].id;
				*os<<setw(columnWidth)<<scores[h][i].label<<setw(columnWidth)<<scores[h][i].score<<"\n";
			}

		}
	}
}

void NTestDistribution::ReportConfig(ostream *os) {
	if (os == NULL)
		os = &cout;
	*os<<"N Distribution Permutation Distribution\n";
	*os<<"Random Seed: "<<seed<<endl;
	*os<<"Test Count: "<<scores[0].size()<<endl;
}

//Assumes that the vector has been sorted
float NTestDistribution::GetPValue(uint modelSize, float score) {
	if (!isSorted) 
		Sort();

	uint i=0;
	uint idx=modelSize-1;
	uint count=scores[idx].size();

	while (scores[idx][i++].score > score && i<count);
	return (float)i/(float)count;
}
void NTestDistribution::Append( float score, uint modelSize, uint test, const char *label) {
	MPILOCK;
	uint idx=modelSize;
	assert(idx < scores.size());
	if (!(test < scores[idx].size()))
		cout<<"!! Trying to add too many test values. test ID="<<test<<" real size="<<scores[idx].size()<<"\n";
	assert(test < scores[idx].size()); 
	isSorted=false;

	if (scores[idx][test].score<score) 
		scores[idx][test].Evaluate(label, score);
		else
			cout<<"Skipping assigning "<<score<<" to "<<test<<" "<<idx<<" to replace "<<scores[idx][test].score<<"\n";
	MPIUNLOCK;

}

void NTestDistribution::Append(float score,  uint test, SnpAligned *snp) {
	MPILOCK;
	if (snp) {
		uint idx=snp->GetLabelCount()-1;
		assert(idx < scores.size());
		if (!(test < scores[idx].size()))
			cout<<"!! Trying to add too many test values. test ID="<<test<<" real size="<<scores[idx].size()<<"\n";
		assert(test < scores[idx].size()); 
		isSorted=false;
		if (scores[idx][test].score<score) 
			scores[idx][test].Evaluate(snp, score);
		else
			cout<<"Skipping assigning "<<score<<" to "<<test<<" "<<idx<<" to replace "<<scores[idx][test].score<<"\n";
	}
	else {
		cout<<"!! An internal error has occurred. The applicated attempted to append a score without an associated model size!\n";
		assert(snp);
	}
	
	MPIUNLOCK;
}

void NTestDistribution::Append(uint test, SnpAligned *snp) {
	MPILOCK;
	if (snp) {
		uint idx=snp->GetLabelCount()-1;
		float score = snp->GetLastMdEval();
		assert(idx < scores.size());
		assert(test < scores[idx].size()); 
		isSorted=false;
		if (scores[idx][test].score<score) 
			scores[idx][test] = snp;
	}
	MPIUNLOCK;
		
}

void NTestDistribution::Sort() {
	for (uint i=0; i<scores.size(); i++) 
		sort(scores[i].begin(), scores[i].end());

	isSorted=true;
}




}
}
