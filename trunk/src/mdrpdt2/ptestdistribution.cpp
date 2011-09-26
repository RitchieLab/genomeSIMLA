//
// C++ Implementation: ptestdistribution
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ptestdistribution.h"
#include <iomanip>

using namespace std;

namespace MdrPDT {

namespace Distribution {

Distribution::Distribution(int testCount) : scores(testCount, ScoreCard("", 0.0)) { }

Distribution::~Distribution() { }

void Distribution::AppendTest(int idx, const char *id, float val) { 
	if (val>scores[idx].score) 
		scores[idx] = ScoreCard(id, val); 
}

float Distribution::Evaluate(float val) { 
	ScoreCard score("", val);
	Init();
	int position = scores.size();
	TDistributionNode *node = distribution.FindNearestMin(score);
	if (node) {
		position = node->GetData();
	}
	return (float)position/(float)scores.size();
}

void Distribution::ReportDistribution(std::ostream& os) {
	Init();
	int pos = 1;
	TDistributionNode *node = distribution.GetFirst();
	os<<setw(8)<<"Index"<<setw(10)<<setprecision(4)<<"MOR"<<setw(10)<<"Model"<<setw(10)<<"R. Index"<<"\n";
	os<<"--------------------------------------\n";
	while (node) {
		os<<setw(8)<<pos++
			<<setw(10)<<setprecision(4)<<node->GetKey().score
			<<setw(10)<<node->GetKey().id
			<<setw(10)<<node->GetData()<<"\n";
		node=node->GetNext();
	}
}

void Distribution::Init() {
	if (distribution.GetCount() == 0) {
		sort(scores.rbegin(), scores.rend());
		
		int count = scores.size();
		for (int i=0; i<count; i++) 
			distribution.Add(scores[i], i+1);
	}
}

int PTestDistribution::GetSignificantDigits() {
	int testCount = GetDistributionSize();
	int sigDig = 0;
	
	while (testCount > 1) {
		sigDig++;
		testCount/=10;
	}
	return sigDig;
}

OmnibusDistribution::OmnibusDistribution(int testCount) : distribution(testCount) { }

void OmnibusDistribution::AppendTest(int modelSize, int idx, const char *id, float val) {
	distribution.AppendTest(idx, id, val);
}

float OmnibusDistribution::Evaluate(int modelSize, float score) {
	return distribution.Evaluate(score);
}

void OmnibusDistribution::ReportDistribution(int modelSize, std::ostream& os) {
	distribution.ReportDistribution(os);
}




}


}
