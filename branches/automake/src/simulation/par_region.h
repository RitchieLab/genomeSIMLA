//
// C++ Interface: par_region
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPAR_REGION_H
#define SIMULATIONPAR_REGION_H

#include "utility/random.h"
#include <math.h>
#include <vector>
#include "locusmanager.h"
//#include "allelesource.h"
#include "utility/rbtree.h"
#include "utility/stat.h"

#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif
namespace Simulation {



/**
	@brief Simple, growable structure used to keep track with one or more PAR region in the Y Chromosomes
	@Usage Basically, ChromPoolXY will feed each consecutive neighboring locus paris in a PAR to the 
		main PAR_Region. It will split into multiple PARs if they exist. 
		To draw a position along the PAR Region, the pool will make a draw against the combined length
		of the PAR(s) and query the region for the converted location (which takes into account the
		physical position of the bounding Loci)

	@author Eric Torstenson
*/
template <class T>
class PAR_Region {
public:

	PAR_Region(LocusManager<T> *loci, typename LocusManager<T>::PositionLookup &yLookup);
	PAR_Region();
	~PAR_Region();
	/**
	 * @brief Inserts a two locus region into the PAR
	 * @note This could grow the current region, or cause a new region to be created (or grown)
	 */
	void InsertRegion(T& start, T& stop);

	/**
	 * @brief Determines the total size of this and child PARs
	 */
	float CalculateSize();
	/**
	 * @brief Converts a local draw to a real location on the chromosome. 
	 * @note this converts a localized draw on the length of all pars into a real position into the chromosome
	 */
//	uint ConvertXO(uint contPoint);

	/**
	 * @brief Returns the phase and sets eventcount to the number of events drawn from this and all subsequent PARs
	 */
//	bool XOEventCount(Utility::Random& rnd, int& eventCount, double pLambda);

	/**
	 * @brief Returns phase for output type and builds the vector of cross-overs within the PAR regions
	 * @note XO is based on the Y genetic positions, since this structure will only be used when crossing over an X and a Y (Male)
	 */
	bool BuildXOEventList(Utility::Random& rnd, std::vector<size_t>& events);
	
	void Init(LocusManager<T> *loci);
protected:
	T start;
	T stop;
	PAR_Region *next;
	LocusManager<T> *loci;
	typename LocusManager<T>::PositionLookup yLookup;				///<Locate SNP based on Y genetic position
};



template<class T>
inline
PAR_Region<T>::PAR_Region(LocusManager<T>* loci, typename LocusManager<T>::PositionLookup& ydata): next(NULL), loci(loci), yLookup(ydata) { }

template<class T>
inline
PAR_Region<T>::PAR_Region() : next(NULL), loci(NULL) { }

template<class T>
inline
void PAR_Region<T>::Init(LocusManager<T> *loci) {
	this->loci = loci;
	yLookup.Clear();

	int count = loci->LocusCount();
	for (int i=0; i<count; i++) {
		T *locus = loci->At(i);
		float position = locus->GetMapPositionY();
		if (position > 0.0000000001) 
			yLookup.Set(position, locus);
	}
}



template<class T>
inline
PAR_Region<T>::~PAR_Region() {
	if (next)
		delete next;
}

template<class T>
inline
void PAR_Region<T>::InsertRegion(T& start, T& stop) {
	//If this is so not part of the local region, we add it to the next
	if (start.GetMapPositionY() > this->stop.GetMapPositionY() && this->stop.GetMapPositionY() != 0) {
		if (next == NULL) 
			next = new PAR_Region(loci, yLookup);
		next->InsertRegion(start, stop);
	} else {
		if (this->start.GetMapPositionY() < this->start.GetMapPositionY() || this->stop.GetMapPositionY() == 0) {
			assert(stop.GetMapPositionY() > this->stop.GetMapPositionY());
			this->start = start;
		}
		if (stop.GetMapPositionY() > this->stop.GetMapPositionY()) {
			this->stop = stop;
		}
	}
		
}


template<class T>
inline
bool PAR_Region<T>::BuildXOEventList(Utility::Random& rnd, std::vector<size_t>& events) {
	//If it dies here, we failed to call Init(locusArray);
static int count =0;

	assert(loci);
	float startPos = start.GetMapPositionY();
	float distance = stop.GetMapPositionY() - start.GetMapPositionY(); 
	float pLambda  = distance * 0.01;

	int eventCount = PoissonEventCount(rnd, pLambda);
	bool phase = eventCount %2;

	for (int i=0; i<eventCount; i++) {
		float loc = rnd(distance) + startPos;
		typename LocusManager<T>::PositionLookupNode *node = yLookup.FindNearestMin(loc);
		T *minLoc = node->GetData();
		//T *minLoc = loci->At(loc);
		int idx = minLoc->GetID();
		vector<size_t>::iterator other = find(events.begin(), events.end(), idx);
		if (other != events.end())
			events.erase(other);
		else
 			events.push_back(idx);
	}
	sort(events.begin(), events.end());
	//This is pretty dumb, but we only care about the first one, so it should be empty to start with
	
	if (next)
		//We can throw the child's phase away, since it's irrelavent
		next->BuildXOEventList(rnd, events);

/*if (count++%100 == 0)
cerr<<"PAR: "<<pLambda<<" "<<phase<<" "<<events.size()<<"\n";
*/
	return phase;
}
/*
template<class T>
inline
bool PAR_Region<T>::XOEventCount(Utility::Random& rnd, int& eventCount, double pLambda) {
	eventCount = PoissonEventCount(rnd, pLambda);
	bool phase = !(eventCount / 2 * 2 == eventCount);

	int childEvents = 0;
	if (next)
		XOEventCount(rnd, childEvents, pLambda);
	eventCount += childEvents;
	return phase;
}

template<class T>
inline
float PAR_Region<T>::CalculateSize() {
	float distance = 0.0;
	if (next) 
		distance = next->CalculateSize();
	distance += stop.GetMapPositionY() * 1000000 - start.GetMapPositionY() * 1000000;
	return distance;
}

template<class T>
inline
uint PAR_Region<T>::ConvertXO(uint xo)  {
	float distance = (stop.GetMapPositionY() - start.GetMapPositionY()) * 1000000;
	if (xo <= distance)
		return start.GetLocation()+xo;
	else if (next) 
		return next->ConvertXO(xo - (distance));

	assert(0);
	return 0;
}*/

#ifdef CPPUNIT
class PARTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( PARTest );
	CPPUNIT_TEST( TestPAR );
	CPPUNIT_TEST_SUITE_END();
public:
	PARTest();
	~PARTest();

	void setUp();
	void tearDown();


	void TestPAR();
};

#endif //CPPUNIT
}

#endif
