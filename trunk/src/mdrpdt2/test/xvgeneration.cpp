//
// C++ Implementation: xvgeneration
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "xvgeneration.h"

namespace MdrPDT {

namespace Test {

bool XvGeneration::TestDistribution() {
	bool success = true;

	MdrPDT::NTestDistribution ntest(1);
	
	std::vector<int> values;
	for (int i=0; i<10; i++) {
		values.push_back(i+1);
	}
	random_shuffle(values.begin(), values.end());
	for (int i=0; i<values.size(); i++) {
		ntest.AppendTest(0, (float)values[i]);
//		ntest.AppendTest(1, (float)(100*values[i]));
//		ntest.AppendTest(2, (float)(values[i]/100));
	}
	ntest.ReportDistribution(0, cout);
	
	float pValue = ntest.Evaluate(0, 4.9);
	cout<<pValue<<" Should be 0.005\n";
	pValue = ntest.Evaluate(0, 5.0);
	cout<<pValue<<" Should be 0.005\n";
	pValue = ntest.Evaluate(0, 0.5);
	cout<<pValue<<" Should be 0.000\n";
	pValue = ntest.Evaluate(0, 1000.0);
	cout<<pValue<<" Should be 1.0\n";
	pValue = ntest.Evaluate(0, 999.55);
	cout<<pValue<<" Should be 0.99\n";
	cout<<(float)(ntest.Evaluate(0, 8.9))<<" Should be 0.2\n";
	cout<<(float)(ntest.Evaluate(0, 9.1))<<" Should be 0.1\n";
	cout<<(float)(ntest.Evaluate(0, 9.9))<<" Should be 0.1\n";
	cout<<(float)(ntest.Evaluate(0, 10.1))<<" Should be 0.1\n";
	cout<<(float)(ntest.Evaluate(0, 0.1))<<" Should be 1.0\n";
	
	return success;
}



}

}
