#ifndef STAT_H
#define STAT_H 
//
// C++ Implementation: stat.h
//
// Description: Various Statistics calculations that have somewhat general applicability
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "utility/random.h"
#include <vector>
#include <math.h>


namespace Utility {

size_t PoissonEventCount(Utility::Random& rnd, double labmda);
/**
 * Various statistics calculations that have somewhat general applicability
 * @param rnd This is the random number generator to be used
 * @param poissonLambda The lambda to be used
 */
inline
size_t PoissonEventCount(Utility::Random& rnd, double lambda) {
	double p=exp(-lambda), 
		g=p, 
		u=rnd.drand();

	size_t k=0;
	while (u>g)		{
		p*=(lambda/(double)(++k));
		g+=p;
	}
	return k;
};

/**
 * @brief Very simple accumulation class for calculating mean and std-deviation
 */
template<class T>
class Accum {
public:
	std::vector<T> data;
	int sum;
	Accum(int size = 0) : sum(0) {
		data.reserve(size);
	}

	void AddValue(T val) {
		data.push_back(val);
		sum+=val;
	}

	float Mean() {
		return sum/(float)data.size();
	}

	T Sum() {
		return sum;
	}

	float StdDev() {
		typename std::vector<T>::iterator itr = data.begin();
		typename std::vector<T>::iterator end = data.end();

		float mean = Mean();

		double vsum = 0.0;
		while (itr != end) {
			double variance = (double)*itr++ - mean;

			vsum += variance*variance;
		}

		return sqrt(vsum/(float)data.size());
	}
};
}

#endif //STAT_H
