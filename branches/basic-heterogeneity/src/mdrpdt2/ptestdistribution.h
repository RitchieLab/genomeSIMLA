//
// C++ Interface: ptestdistribution
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTPTESTDISTRIBUTION_H
#define MDRPDTPTESTDISTRIBUTION_H
#include "utility/types.h"
#include "utility/rbtree.h"
#include <string>

namespace MdrPDT {

namespace Distribution {

struct ScoreCard {
	std::string id;
	float score;
	ScoreCard() : id(""), score(0.0) { }
	ScoreCard(std::string id, float score) : id(id), score(score) { }
	bool operator<(const ScoreCard& other) const { return score<other.score; }
};

struct ScoreEval {
	int operator()(const ScoreCard& l, const ScoreCard& r) const {
		if (l.score<r.score) return -1;
		if (l.score>r.score) return 1;
		return 0;
	}
};



typedef Utility::RBTree<ScoreCard, int, ScoreEval> TDistribution;
typedef Utility::RBTreeNode<ScoreCard, int, ScoreEval> TDistributionNode;


/**
@Brief Basic Distribution functionality

	@author 
*/

/**
 * @brief This keeps track of a single distribution..provides basic lookup functionality
 */
class Distribution{
public:

    Distribution(int testCount);
    ~Distribution();
	void AppendTest(int idx, const char *id, float val);
	
	float Evaluate(float val);
	void ReportDistribution(std::ostream& os);


	int GetDistributionSize() { return distribution.GetCount(); }
protected:
	void Init();

	std::vector<ScoreCard> scores;
	TDistribution distribution;

};

/**
 * @brief Base class for the distrbution classes
 */
class PTestDistribution {
public:
	virtual void AppendTest(int modelSize, int idx, const char *id, float val)=0;
	virtual float Evaluate(int modelSize, float score)=0;
	virtual void ReportDistribution(int modelSize, std::ostream& os)=0;
	virtual int GetDistributionSize()=0;
	int GetSignificantDigits();
	virtual ~PTestDistribution() { }
};

class OmnibusDistribution: public PTestDistribution {
public:
	OmnibusDistribution(int testCount);

	/**
	 * @brief Append a test to the distribution
	 * @param modelSize The order of the model being added (doesn't matter for omnibus, but...part of the generic distribution model
	 * @param idx The position of the model (this allows us to keep only the best for omnibus)
	 * @param id The model ID...useful for logging
	 * @param val The actual score being stored
	 */
	void AppendTest(int modelSize, int idx, const char *id, float val);

	/**
	 * @brief Evaluate a given score for it's p-Value
	 */
	float Evaluate(int modelSize, float score);

	/**
	 * @brief Write the distribution to the stream (for a given order)
	 */
	void ReportDistribution(int modelSize, std::ostream& os);

	int GetDistributionSize() { return distribution.GetDistributionSize(); }
protected:
	Distribution distribution;				///<The actual distribution
};

class NTestDistribution: public PTestDistribution {
public:
	NTestDistribution(int modelCount, int testCount) : distributions(modelCount, testCount) { }
	
	/**
	 * @brief Append a test to the distribution
	 * @param modelSize The order of the model being added 
	 * @param idx The position of the model (this allows us to keep only the best for omnibus)
	 * @param id The model ID...useful for logging
	 * @param val The actual score being stored
	 */	
	void AppendTest(int modelSize, int idx, const char *id, float val) {
		distributions[modelSize].AppendTest(idx, id, val);
	}
	/**
	 * @brief Evaluate a given score for it's p-Value
	 */
	float Evaluate(int modelSize, float value) {
		return distributions[modelSize].Evaluate(value); 
	}

	/**
	 * @brief Write the distribution to the stream (for a given order)
	 */
	void ReportDistribution(int modelSize, std::ostream& os) {
		distributions[modelSize].ReportDistribution(os);
	}
		int GetDistributionSize() { return distributions[0].GetDistributionSize(); }
protected:
	std::vector<Distribution> distributions;

};
}
}

#endif
