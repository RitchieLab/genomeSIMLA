//
// C++ Implementation: blocklistnode
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iomanip>
#include "blocklistnode.h"

//#define DEBUG_LD 1
#ifdef DEBUG_LD
#define LD_WRITE(A) cout<<A
#else
#define LD_WRITE(A) 
#endif

#include <sstream>

namespace Simulation {

BlockListNode::BlockListNode(float criterion) : previous(NULL), next(NULL), completedList(NULL), deletedList(NULL), firstIdx(0), lastIdx(0), criterion(criterion), avgMAF(0), totalDistance(0), 
		blockDensity(0), minAlleleFreq(NULL), maxAlleleFreq(NULL), previousEdge(0) { }

BlockListNode::BlockListNode(float criterion, BlockListNode *previous, BlockList *completedList, BlockList *deletedList) : previous(previous), next(NULL), completedList(completedList), deletedList(deletedList), firstIdx(0), lastIdx(0), criterion(criterion), avgMAF(0), totalDistance(0), 
		blockDensity(0), minAlleleFreq(NULL), maxAlleleFreq(NULL), previousEdge(0) { }

BlockListNode::BlockListNode(const BlockListNode& other) :	loci(other.loci), previous(other.previous),
		 next(other.next), completedList(other.completedList), deletedList(other.deletedList), firstIdx(other.firstIdx), lastIdx(other.lastIdx), criterion(other.criterion), avgMAF(other.avgMAF), totalDistance(other.totalDistance), 
		blockDensity(other.blockDensity), minAlleleFreq(other.minAlleleFreq), maxAlleleFreq(other.maxAlleleFreq), previousEdge(other.previousEdge) { }
	

void BlockListNode::SetIndexBounds(uint firstIdx, uint lastIdx) {
	this->firstIdx = firstIdx;
	this->lastIdx  = lastIdx;
}

BlockListNode::~BlockListNode()
{

	
//	cout<<"Bye Bye\n";
}

size_t BlockListNode::GetFirstIndex() const {
	return firstIdx;
}

size_t BlockListNode::GetLastIndex() const {
	return lastIdx;
}

size_t BlockListNode::GetBlockSize() const {
	return loci.size();
}

/**
 * @brief Returns the block density
 */
size_t BlockListNode::GetBlockDensity() const {
	return (size_t)blockDensity;
}




void BlockListNode::CompleteNode(vector<Locus> &loci) {
#ifdef DEBUG_LD
	const char *labelA = loci[firstIdx].GetLabel().c_str();
	const char *labelB = loci[lastIdx].GetLabel().c_str();
#endif
	if (previous)
		previous->next = next;
	if (next) {		
		next->previous = previous;
	}
	
	LD_WRITE("\n >  "<<setw(10)<<firstIdx<<" "<<labelA<<" : "<<setw(10)<<lastIdx<<" "<<labelB<<" ("<<previousEdge<<") ");
	//Just some quick sanity check. If the bounds are senseless, let's trash it
	if (firstIdx < lastIdx && lastIdx > previousEdge) 	{
		LD_WRITE(setw(10)<<"Kept!"<<" "<<completedList->size()<<"\n");
		blockDensity = (loci[lastIdx].GetLocation() - loci[firstIdx].GetLocation() ) / (lastIdx - firstIdx);

		if (next)
			next->previousEdge = lastIdx;

assert(lastIdx - firstIdx < 1000);
		for (uint i = firstIdx; i<=lastIdx; i++) 
			this->loci.push_back(loci[i]);

		completedList->push_back(this);
	}
	else {
		lastIdx = firstIdx = 0;
		LD_WRITE(setw(10)<<"Deleted\n");
		if (next) {
			if (next->previousEdge < previousEdge) {
//				cout<<" Setting previous Edge for "<<next->firstIdx<<"x"<<next->lastIdx<<"\n";
				next->previousEdge = previousEdge;
			}
		}

		deletedList->push_back(this);
	}
	next = previous = NULL;

}

string BlockListNode::GetLabel() {
	const char *labelA = loci[0].GetLabel().c_str();
	const char *labelB = loci[loci.size() - 1].GetLabel().c_str();

	stringstream ss;
	ss<<labelA<<"_"<<labelB;
	return ss.str();
}

void BlockListNode::DeleteNode(vector<Locus> &loci) {
#ifdef DEBUG_LD
	const char *labelA = loci[firstIdx].GetLabel().c_str();
	const char *labelB = loci[lastIdx].GetLabel().c_str();
#endif
	if (previous)
		previous->next = next;
	if (next)
		next->previous = previous;
	LD_WRITE("\n -> "<<setw(10)<<firstIdx<<" "<<labelA<<" : "<<setw(10)<<lastIdx<<" "<<labelB<<" ("<<previousEdge<<") Deleted");
	if (next) {
		if (next->previousEdge < previousEdge) {
//			cout<<" Setting previous Edge for "<<next->firstIdx<<"x"<<next->lastIdx<<"\n";
			next->previousEdge = previousEdge;
		}
	}
	deletedList->push_back(this);
	next = previous = NULL;
}

/**
 * Various situations:
 *    	- We reach a condition where the current node is recognizably finished
 *    		- idxA > lastIdx
 *          - idxA == lastIdx - 1 && idxB == lastIdx
 *          - idxA + 1 ==  lastIdx == idxB - 2
 *    	- We reach a value that meets our criterion- we can just move on to the next node
 * 		- We reach a value that fails to meet our criterion, so we need to do one of the following:
 * 			- if (idxA == firstIndex && idxB == firstIndex + 1)  -- We can get rid of the left most index since it's root is junk
 * 			- else if (idxA >= firstIndex ) 					 -- We either perform next->firstIdx++ or we create a new node at next (idxA + 1, lastIdx)
 */
void BlockListNode::TruncateBlock(uint idxA, uint idxB) {
	if (idxB <= lastIdx) {
		if (next == NULL)
			if (idxA + 1< lastIdx)  {
				BlockListNode *newNode = NewNode(this, idxA + 2, lastIdx);
				next = newNode;
			}
		lastIdx = idxB - 1;
	}
}
void BlockListNode::Append(uint idxA, uint idxB, vector<Locus> &loci, float value) {
	//If this is true, the LD values are for peripheral areas we aren't interested in
	//cout<<"-"<<firstIdx<<" - "<<lastIdx<<" ("<<idxA<<"x"<<idxB<<") = "<<value<<"\n";
	BlockListNode *nextSet = next;
#ifdef DEBUG_LD
	const char *labelA = loci[idxA].GetLabel().c_str();
	const char *labelB = loci[idxB].GetLabel().c_str();
	const char *blStart = loci[firstIdx].GetLabel().c_str();
	const char *blStop  = loci[lastIdx].GetLabel().c_str();


#endif

	if (firstIdx == lastIdx) {
		DeleteNode(loci);
		return;
	}

	if (idxA < firstIdx) 
		return;
	
	LD_WRITE("Considering LD Unit: "<<idxA<<"x"<<idxB<<" "<<labelA<<" x "<<labelB<<" ("<<blStart<<" x "<<blStop<<") "<<value<<" : ");
	
	//If we are sitting on the last node of potential block, 
	//evaluate it and send the block for completion	
	if (idxA == lastIdx - 1) {
		//assert(next != NULL);
		if (value < criterion) {
			LD_WRITE(" - Truncated\n");
			lastIdx--;
			CompleteNode(loci);
		}
		else {
			LD_WRITE(" - Remains intact\n");
			CompleteNode(loci);
		}
	}
	else if (idxA >= lastIdx) {
		LD_WRITE(" - Exceeded whole block size\n");
		//cout<<"IdxA >= lastIdx ("<<idxA<<" >= "<<lastIdx<<") Sending for completion\n";
		CompleteNode(loci);
	}
	else if (idxA >= firstIdx && idxB <= lastIdx) {
		LD_WRITE(" - Within Bounds ");
		//This is the last corner of the region
		if (((idxA + 1) == lastIdx) && (lastIdx == idxB)) {
			LD_WRITE(" RH Position");
			//cout<<"This is the last value in the current block: ("<<firstIdx<<" "<<lastIdx<<") "<<value<<" : "<<idxA<<" "<<idxB<<"\n";
			if (value >= criterion)  {
				LD_WRITE(" Accepted\n");
				//cout<<"Completing the node\n";
				CompleteNode(loci);	
			}
			else {
				LD_WRITE(" Truncated\n");
				//cout<<"Trashing it \n";
				lastIdx--;
				CompleteNode(loci);
			}
		} else {	
			if (value < criterion) {
				LD_WRITE(" Middle. Truncation required ("<<setprecision(4)<<value<<"<"<<setprecision(4)<<criterion<<") ");
				//Check to see if the previous node(s) were good
				if (firstIdx < idxA) {
					LD_WRITE(" Keeping previous ");
					if (next == NULL)
						if (idxA + 1< lastIdx)  {
	 						BlockListNode *newNode = NewNode(this, idxA + 1, lastIdx);
							nextSet = next = newNode;
							//cout<<"New Node: "<<next->firstIdx<<" x "<<next->lastIdx<<"\n";
						}
					else
						next->firstIdx = idxA + 1;
					
					//In this case, the previous node(s) were good.
					lastIdx = idxB - 1;
					LD_WRITE("-> ("<<loci[firstIdx].GetLabel()<<" x "<<loci[lastIdx].GetLabel()<<")");
					if (idxA >= lastIdx)
						CompleteNode(loci);
				}	
				else {				
					LD_WRITE(" Discarding previous\n");
					//corner...which is junk...
					if (idxA == idxB -1) {
						firstIdx++;
						if (next)
							next->firstIdx++;
					}
					//Let's perform a split and check for completion
					else {
						LD_WRITE(" Splitting ");
						if (next == NULL && idxA + 1< lastIdx)  {
							BlockListNode *newNode = NewNode(this, idxA + 1, lastIdx);
							nextSet = next = newNode;
							//cout<<"New Node: "<<next->firstIdx<<" x "<<next->lastIdx<<"\n";
						}	
						next->firstIdx = idxA + 1;
						lastIdx = idxB - 1;
						LD_WRITE("-> ("<<loci[firstIdx].GetLabel()<<" x "<<loci[lastIdx].GetLabel()<<")");
						//If we reach this point, the end of the triangle has already been proven
						if (idxA >= idxB - 2) 
							CompleteNode(loci);
					}
				}
			}
		}
	}
	LD_WRITE("\n");
	if (nextSet)
		nextSet->Append(idxA, idxB, loci, value);	
}

/**
 * @brief returns the number of snps contained within the block
 */
size_t BlockListNode::BlockCount() const {
	return loci.size();
}

/**
 * @brief REturns the average minor allele frequency for the block
 */
double BlockListNode::GetAverageMAF() const {
	return avgMAF;
}

/**
 * @brief Return the first location on the block
 */
uint BlockListNode::GetFirstLocation() const {
	return firstLocation;
}

/**
 * @brief returns the last location on the bloc, 
 */
uint BlockListNode::GetLastLocation() const {
	return lastLocation;
}
/*
bool BlockListNode::operator<(const BlockListNode& other) const {
	if (selectionScore == other.selectionScore)
		return loci.size() < other.loci.size();
	else
		return selectionScore < other.selectionScore;
}*/


}
