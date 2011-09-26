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

namespace Genetics {

namespace Reporting {


inline
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
			*os<<setw(columnWidth)<<"Test"<<setw(columnWidth)<<"Model"<<setw(columnWidth)<<""<<endl;
			*os<<setw(columnWidth)<<"Number"<<setw(columnWidth)<<" ID "<<setw(columnWidth)<<"Score"<<endl;
			for (uint i=0; i<count; i++) { 
				if (scores[h][i].snp) {
					*os<<setiosflags(ios::fixed | ios::showpoint) <<setprecision(3);
					*os<<setw(columnWidth)<<i<<setw(columnWidth)<<scores[h][i].snp->GetLabel()<<setw(columnWidth)<<scores[h][i].score<<"\n";
				}
			}
		}
	}
}

inline
void NTestDistribution::ReportConfig(ostream *os) {
	if (os == NULL)
		os = &cout;
	*os<<"N Distribution Permutation Distribution\n";
	*os<<"Random Seed: "<<seed<<endl;
	*os<<"Test Count: "<<scores[0].size()<<endl;
}

//Assumes that the vector has been sorted
inline
float NTestDistribution::GetPValue(uint modelSize, float score) {
	if (!isSorted) 
		Sort();

	uint idx=modelSize-1;

	uint count=scores[idx].size();
	float pValue=1.0;

	for (uint i=count; i>0; i--){ 
		if (scores[idx][i-1].score<score) {		
			pValue=(float)(count-i+1)/(float)count;
			break;
		}
	}
	return pValue;
}

inline
void NTestDistribution::Append(float score,  uint test, SnpAligned *snp) {
	if (snp) {
		uint idx=snp->GetLabelCount()-1;
		assert(idx < scores.size());
		assert(test < scores[idx].size()); 
		isSorted=false;
		if (scores[idx][test].score<score) 
			scores[idx][test].Evaluate(snp, score);
	}
}


void NTestDistribution::Append(uint test, SnpAligned *snp) {
	if (snp) {
		uint idx=snp->GetLabelCount()-1;
		float score = snp->GetLastMdEval();
		isSorted=false;
		if (scores[idx][test].score<score) 
			scores[idx][test] = snp;
	}
		
}

void NTestDistribution::Sort() {
	for (uint i=0; i<scores.size(); i++) 
		sort(scores[i].begin(), scores[i].end());

	isSorted=true;
}

}
}
