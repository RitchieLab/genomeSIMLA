//
// C++ Interface: ptestdistribution
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_REPORTINGPTESTDISTRIBUTION_H
#define GENETICS_REPORTINGPTESTDISTRIBUTION_H

#include <string>
#include <vector>
#include "snpaligned.h"
#include <math.h>
namespace Genetics {

namespace Distribution {


using namespace std;

/**
 * Contents of the distribution
 */
struct DistContents {
	string modelLabel;				///<Label of the model with this score. Just for dumping distribution
	float value;					///<The value at a given test
	int testID;						///<Used to keep up with which test a given entry is associated with

	DistContents(const char *lbl, float val, int testID) : modelLabel(lbl), value(val), testID(testID) {}
	DistContents() : modelLabel("nomodel"), value(0.0), testID(0) { }

	/**
	 * We are going to sort this according to largest value
	 */
	bool operator<(const DistContents& other) const { return value > other.value; }
};

/**
 * Base class for the distributions
 */
class PTestDistribution {
public:
	PTestDistribution(const char *label, int seed=1313) : label(label), isSorted(false), seed(seed) {}
	PTestDistribution(int seed=1313) : label(""), isSorted(false), seed(seed) {}
	virtual ~PTestDistribution() {}
	/**
	 * Explicit assignment - Only valid entries are to be appended
	 */
	virtual void AppendTest(float value, uint modelSize, const char *label, int testID=-1) = 0;

	/**
	 * Assignment - The other tidbits are found in the snp
	 */
	virtual void AppendTest(float value, SnpAligned *snp, int testID=-1) = 0;

	/**
	 * Return the pvalue associated with value
	 */
	virtual float GetPValue(float value, uint modelSize) = 0;

	virtual string GetLabel();

	/**
	 * Sorts the entries
	 */
	virtual void Sort() = 0;

	/**
	 * Write the distribution to a stream
	 */
	virtual void DumpDistribution(ostream* os) = 0;

	/**
	 * Read the distribution from a file
	 */
	virtual void LoadDistribution(istream *os) = 0;

	/**
	 * Returns the number of tests associated with this distribution
	 */
	virtual uint GetTestCount(uint modelSize)=0;

	/**
	 * Since we might have more than one type of distribution, we will label the field
	 */
	void SetValueLabel(const char *label) { 
		this->label = label;
	}
	
protected:
	string label;								///<Used to describe the type of distribution
	bool isSorted;								///<Used to determine if the test is sorted
	int seed;									///<Record the seed being used
};

/**
NTest Distribution - This distribution maintains separate distributions for each order of models considered

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class NTest : public PTestDistribution {
public:
    NTest(uint maxModelSize, int seed, char *label);
    virtual ~NTest();

	void AppendTest(float value, uint modelSize, const char *modelLabel, int testID=-1);
	void AppendTest(float value, SnpAligned *snp, int testID=-1);
	float GetPValue(float value, uint modelSize);
	uint GetTestCount(uint modelSize);
	void DumpDistribution(ostream* os);
	void LoadDistribution(istream *os);
	void Sort();
protected:
	/**
	 * Values. distribution[order][testnumber] = score
	 */
	vector<DistContents> *distribution;	

	/**
	 * Keep up with the size of the array of distributions
	 */
	uint maxModelSize;
};

/**
 * Omnibus Distribution - This keeps a single distribution regardless of the order
 */
class Omnibus : public PTestDistribution {
public:
	Omnibus(int testCount, int seed, const char *label, float defInitValue = 0.0);
	~Omnibus();
	void AppendTest(float value, uint modelSize, const char *modelLabel, int testID=-1);
	void AppendTest(float value, SnpAligned *snp, int testID=-1);
	float GetPValue(float value, uint modelSize);
	uint GetTestCount(uint modelSize);
	void Sort();
	void DumpDistribution(ostream* os); 
	void LoadDistribution(istream *os);
protected:

	/**
	 * Values. distribution[testnumber] = score
	 */
	vector<DistContents> distribution;	
	int lastID;

	int testCount;
};

class InvertedOmnibus : public Omnibus {
public:
	InvertedOmnibus(int testCount, int seed, const char *label, float initValue) : Omnibus(testCount, seed, label, initValue) {	}
	void AppendTest(float value, uint modelSize, const char *modelLabel, int testID=-1);
	void AppendTest(float value, SnpAligned *snp, int testID=-1);
	void Sort();	
	float GetPValue(float value, uint modelSize);
};

inline
string PTestDistribution::GetLabel() {
	return label;
}

}

}

#endif
