//
// C++ Interface: kinshipcalculator
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONKINSHIPCALCULATOR_H
#define SIMULATIONKINSHIPCALCULATOR_H
#include <set>
#include <vector>
#include <string>
#include <map>
#include "utility/types.h"
#include "utility/rbtree.h"
#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif

namespace Simulation {
using namespace std;
/**
@brief Calculates the kinship for a pair of chromosomes based on the allele source tree

	@author Eric Torstenson
*/
template <class T>
class KinshipCalculator{
friend class KinshipCalculatorTest;
public:

typedef Utility::RBTree<size_t, T> SourceTree;
typedef Utility::RBTreeNode<size_t, T> SourceNode;

	class Segment {
	public:
		int start;					///<This represents the start of the inclusive range
		int stop;					///<The end of the inclusive range
		
		/**
	 	 * @brief Returns the amount of overlap, and indicates the direction to move in (L/-1, r/1)
	 	 */
		Segment(int start, int stop) : start(start), stop(stop) { }
		int EvaluateOverlap(const Segment& other, int& direction) {
			if (stop > other.stop)
				direction = 1;
			else
				direction = -1;
			if (start > other.stop || stop < other.start) 
				return 0;

			int segStart = start;
			if (other.start > start)
				segStart = other.start;
			int segEnd = stop;
			if (other.stop < stop)
				segEnd = other.stop;
			return segEnd-segStart+1;
		}

		bool operator==(const Segment& other) {
			return start==other.start && stop==other.stop;
		}
	/**
	 * @brief streaming operator
	 */
    friend std::ostream & operator << (std::ostream & os, const KinshipCalculator::Segment& seg) {
		os<<"Segment("<<seg.start<<","<<seg.stop<<")\n";
		return os;
	}

	};


	typedef vector<typename KinshipCalculator<T>::Segment> SegmentVector;
	typedef map<T, SegmentVector > SegmentMap;


    KinshipCalculator() { }
    ~KinshipCalculator() { }

	/**
	 * @brief Returns the number of overlapping alleles
	 */
	int EvaluateKinship(SourceTree& left, SourceTree& right);

protected:
	/**
	 * @brief Compare two vectors of Segments (assuming that they are both associated with the same key)
	 */
	int Evaluate(SegmentVector& left, SegmentVector& right);
	int Evaluate(SegmentMap& left, SegmentMap& right);
	/**
	 * @brief Used to create the segment map
	 */
	void SegmentTree(SourceTree& tree, SegmentMap& segments);
	void GetKeys(SegmentMap& data, set<T>& keys);
	
};





using namespace Utility;
template <class T>
void KinshipCalculator<T>::SegmentTree(SourceTree& source, SegmentMap& segments) {
	SourceNode* itr = source.GetFirst();
	segments.clear();
	if (itr) {
		int curIndex = itr->GetKey();
		T curSource = itr->GetData();
		itr = itr->GetNext();
		while (itr) {
			T nextSource = itr->GetData();
			int nextIndex = itr->GetKey();
			segments[curSource].push_back(Segment(curIndex, nextIndex-1));
			curIndex = nextIndex;
			itr = itr->GetNext();
			curSource = nextSource;
		}
	}
		
}

template <class T>
void KinshipCalculator<T>::GetKeys(SegmentMap& data, set<T>& keys) {
	typename SegmentMap::iterator itr = data.begin();
	typename SegmentMap::iterator end = data.end();
	keys.clear();
	while (itr != end) {
		keys.insert(itr->first);
		itr++;
	}
}

template <class T>
int KinshipCalculator<T>::EvaluateKinship(SourceTree& left, SourceTree& right) {
	SegmentMap leftMap;
	SegmentMap rightMap;

	SegmentTree(left, leftMap);
	SegmentTree(right, rightMap);
	
	return Evaluate(leftMap, rightMap);
}

template <class T>
int KinshipCalculator<T>::Evaluate(SegmentMap& left, SegmentMap& right) {
	int overlap = 0;

	set<T> leftKeys;
	set<T> rightKeys;
	set<T> intersection;
	
	GetKeys(left, leftKeys);
	GetKeys(right, rightKeys);
	
	set_intersection(leftKeys.begin(), leftKeys.end(), rightKeys.begin(), rightKeys.end(), inserter(intersection, intersection.begin()));

	//OK, now we know which ones to compare:
	typename set<T>::iterator itr = intersection.begin();
	typename set<T>::iterator end = intersection.end();

	while (itr != end) {
		overlap+=Evaluate(left[*itr], right[*itr]);
		itr++;
	}
	return overlap;
}

template <class T>
int KinshipCalculator<T>::Evaluate(SegmentVector& left, SegmentVector& right) {	
	typename vector<typename KinshipCalculator<T>::Segment>::iterator lItr = left.begin();
	typename vector<typename KinshipCalculator<T>::Segment>::iterator lEnd = left.end();
	typename vector<typename KinshipCalculator<T>::Segment>::iterator rItr = right.begin();
	typename vector<typename KinshipCalculator<T>::Segment>::iterator rEnd = right.end();
	
	int direction; 			///< -1 means left needs to advance, 1 means right
	int totalOverlap = 0;
	do {
		totalOverlap += lItr->EvaluateOverlap(*rItr, direction);
		if (direction < 0)
			lItr++;
		else {
			assert(direction > 0);
			rItr++;
		}
	} while (lItr != lEnd && rItr != rEnd);
	return totalOverlap;
}

#ifdef CPPUNIT
class KinshipCalculatorTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( KinshipCalculatorTest );
	CPPUNIT_TEST( TestSegment );
	CPPUNIT_TEST( TestSegmentation );
	CPPUNIT_TEST( TestKinship );
	CPPUNIT_TEST_SUITE_END();
public:
    KinshipCalculatorTest();
    ~KinshipCalculatorTest();


	void setUp();
	void tearDown();
	
	void TestSegment();
	void TestKinship();
	void TestSegmentation();
protected:
	KinshipCalculator<string>::SourceTree left;
	KinshipCalculator<string>::SourceTree right;
};
#endif

}


#endif
