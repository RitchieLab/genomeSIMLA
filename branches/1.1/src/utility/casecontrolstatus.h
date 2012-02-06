//
// C++ Interface: casecontrolstatus
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UTILITYCASECONTROLSTATUS_H
#define UTILITYCASECONTROLSTATUS_H

#include "types.h"
#include <vector>
#include <iostream>
#include "random.h"


namespace Utility {

using namespace std;

typedef size_t size_type;


/**
@brief Stores the affected and unaffected status as two bitset arrays and exposes a function produce randomly generated clones (sets that have the same number of affected and unaffected, but randomly arranged.)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
struct CaseControlStatus{
	BitSetType affected;				///<Mask for affected individuals
	BitSetType unaffected;				///<Mask for unaffected individuals
	BitSetType total;					///<Mask for the all individuals
	size_t affectedWeight;				///<Used to keep track of the weights to be used in the calculations
	size_t unaffectedWeight;			///<Used to keep track of the weights to be used in the calculations
	size_t affectedCount;				///<Cache the counts in case we just need them
	size_t unaffectedCount;				///<Cache the counts in case we just need them
	size_t totalCount;					///<Cache the counts 


	CaseControlStatus(BitSetType &aff) : affected(aff), affectedWeight(1), unaffectedWeight(1) {
		unaffected = ~affected;
		total=affected|unaffected;
		PerformCounts();
	}

	CaseControlStatus(uint indCount =0): affectedWeight(1), unaffectedWeight(1) {
		affected.resize(indCount, false);
		unaffected.resize(indCount, false);
		total.resize(indCount, false);
		PerformCounts();
	}

	CaseControlStatus(BitSetType &aff, BitSetType &unaff, BitSetType &total) : affected(aff), unaffected(unaff), total(total), affectedWeight(1), unaffectedWeight(1) {
		PerformCounts();
	}

	CaseControlStatus(BitSetType &aff, BitSetType &unaff) : affected(aff), unaffected(unaff), affectedWeight(1), unaffectedWeight(1) {
		total=affected | unaffected;
		PerformCounts();
	}
    ~CaseControlStatus() {}
	
	void PerformCounts() {
		//cout<<"Performing counts\n";
		affectedCount = affected.count();
		unaffectedCount = unaffected.count();
		totalCount = total.count();	
	}

	void Reset() {
		affected.clear();
		unaffected.clear();
		total.clear();
		PerformCounts();
	}

	void Flip() {
		affected = ~affected;
		unaffected = ~unaffected;
		total = ~total;
	}

	void Resize(uint newSize) {
		affected.resize(newSize);
		unaffected.resize(newSize);
		total.resize(newSize);
	}

	/**
	 * @brief Creates a copy with randomly distributed affected/unaffected
	 * The distribution matches affected and unaffected counts
	 */
	CaseControlStatus MakeRandomCopy() {
		//Create the bitsets to recieve shuffled content
		BitSetType aff(total.size());
		BitSetType unaff(total.size());
	
		vector<size_type> positions;			//This will hold the positions
		size_type pos = total.find_first();		//The "iterator"
		
		//Iterate over the bitset, recording the positions where we have found a true
		while (pos != BitSetType::npos) {
			positions.push_back(pos);
			pos=total.find_next(pos);
		}
		//Shuffle
		random_shuffle(positions.begin(), positions.end(), Utility::Random::globalGenerator);
	
		size_type affectedCount = affected.count(); 
		//Step through the shuffled vector and set the positions found within
		for (size_type i=0; i<affectedCount; i++)
			aff[positions[i]] = true;
		
		//Let's aquire the unaffected mask
		unaff = total & ~aff;
		CaseControlStatus newbie(aff, unaff, total);
		//return the new item
		return newbie;
	}

	void SetStatus(uint idx, bool isAffected) {
		if (idx >= affected.size() || idx >=unaffected.size()) {
			uint newSize=idx+1;
			affected.resize(newSize);
			unaffected.resize(newSize);
			total.resize(newSize);
		}
		if (isAffected) { 
			if (!affected[idx]) {
				affected[idx]=true;
				affectedCount++;
				totalCount++;
			}
		}
		else {
			if (!unaffected[idx]) {
				unaffected[idx]=true;
				unaffectedCount++;
				totalCount++;
			}
		}
		total[idx]=true;
	}
	void ForceAppend(CaseControlStatus &other) {
		if (total.count() < 1) {
			affected = other.affected;
			unaffected = other.unaffected;
			total=other.total;
			affectedWeight = other.affectedWeight;
			unaffectedWeight = other.unaffectedWeight;
		}
		else {
			if (affected.size() > other.affected.size()) {
				other.affected.resize(affected.size()); 
				other.unaffected.resize(unaffected.size());
				other.total.resize(total.size());
			} else if (affected.size() < other.affected.size()) {
				affected.resize(other.affected.size());
				unaffected.resize(other.unaffected.size());
				total.resize(other.total.size());
			}
				
			affected = affected | other.affected;
			unaffected = unaffected | other.unaffected;
			total = affected | unaffected;
			PerformCounts();
		}
	}		

	void AppendStatus(CaseControlStatus& other) {
		if (total.count() < 1) {
			affected = other.affected;
			unaffected = other.unaffected;
			total=other.total;
			affectedWeight = other.affectedWeight;
			unaffectedWeight = other.unaffectedWeight;
		}
		else {
				
			assert(affectedWeight==other.affectedWeight && unaffectedWeight==other.unaffectedWeight);
			if (affected.size() > other.affected.size()) {
				other.affected.resize(affected.size()); 
				other.unaffected.resize(unaffected.size());
				other.total.resize(total.size());
			} else if (affected.size() < other.affected.size()) {
				affected.resize(other.affected.size());
				unaffected.resize(other.unaffected.size());
				total.resize(other.total.size());
			}
				
			affected = affected | other.affected;
			unaffected = unaffected | other.unaffected;
			total = affected | unaffected;
		}
	}
};


typedef vector<CaseControlStatus> StatusMap;
struct StatusContainer {
	StatusMap status;
	typedef StatusMap::iterator Iterator;
	
	/**
	 * @brief Returns a status that represents all flattened into a single status object.
	 * This is useful for situations where weighted status evaluations aren't required
	 */
	CaseControlStatus CombinedStatus();

	void Resize(uint newSize) {
		size_t stCount=status.size();

		//Let's grab the first one and start there
		for (uint i=0; i<stCount; i++) {
			status[i].affected.resize(newSize);
			status[i].unaffected.resize(newSize);
			status[i].total.resize(newSize);
		}
	}

	void AppendStatus(CaseControlStatus& stat) {
		status.push_back(stat);
	}
		
	void GenerateReport(ostream *os);

	int GetStatusCount(BitSetType &data, bool isAffected);
		
};


inline
CaseControlStatus StatusContainer::CombinedStatus() {
	CaseControlStatus statusTotal;
	size_t count=status.size();

	//Let's grab the first one and start there
	for (uint i=0; i<count; i++) {
		statusTotal.ForceAppend(status[i]);
	}
	return statusTotal;
}

inline
int StatusContainer::GetStatusCount(BitSetType &data, bool isAffected) {
	BitSetType matches; 
	int count=0;
	size_t stCount=status.size();

	for (uint i=0; i<stCount; i++) {
		if (isAffected)
			matches=status[i].affected & data;
		else
			matches=status[i].unaffected & data;
		count+=matches.count();
	}
	return count;
}	

}

#endif
