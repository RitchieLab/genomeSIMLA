//
// C++ Interface: blocklistnode
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONBLOCKLISTNODE_H
#define SIMULATIONBLOCKLISTNODE_H
#include <vector>
#include <map>
#include "locus.h"

namespace Simulation {

using namespace std;


class BlockFilter {
public:
	BlockFilter(bool isValid, int score) : isValid(isValid), score(score) { }	
	virtual ~BlockFilter() {}
	virtual bool IsValid() = 0;
	virtual int GetScore() = 0;
	virtual void Reset() { isValid = false; score = 0; }
	virtual bool Evaluate(Locus &l) = 0;

protected:
	bool isValid;
	int score;

};


class BlockListNode;
typedef vector<BlockListNode *> BlockList;




/**
@brief Used to manage haplotype blocks 
@NOTE Loci order should be based on location such that the snp at IDX 3 will exist to the left of  IDX 5

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BlockListNode{
public:
    BlockListNode(float criterion, BlockListNode *previous, BlockList *completedList, BlockList *deletedList);
	BlockListNode(const BlockListNode& other);
    virtual ~BlockListNode();

	/**
	 * @brief Attempt to assign an LD value to the block
	 * @note If the LD value doesn't fit the necessary criterion, it might cause the current block to spilt
	 */
	virtual void Append(uint idxA, uint idxB, vector<Locus> &loci, float value);
	virtual void TruncateBlock(uint idxA, uint idxB);
	virtual bool ValidateBlocks(vector<Locus> &loci, BlockFilter *filter)=0;

	void SetIndexBounds(uint firstIdx, uint lastIdx);
	
	BlockListNode& operator=(const BlockListNode&other) { assert(0); }
	
	/**
	 * @brief Return a copy of this object without knowing it's exact class
	 */
	virtual BlockListNode *Clone() = 0;
	/**
	 * @brief Returns the block density
	 */
	size_t GetBlockDensity() const;

	size_t GetBlockSize() const;

	/**
	 * @Brief starts the summary (in HTML format)
	 */
	virtual void HtmlSummaryBegin(ostream &os) =0;	

	/**
	 * @brief write block details to the stream (HTML format)
	 * @param os Stream to be written to
	 * @param linkToDetails The name of the detailed plot to be imbedded in the report
	 */
	virtual void HtmlSummary(ostream &os, const char *linkToDetails) =0;
		
	/**
	 * @brief Close out the block's report
	 */
	virtual void HtmlSummaryEnd(ostream &os) =0;

	/**
	 * @brief Sets up the detailed hmtl report
	 * @return the name of the anchor to which the page can be linked
	 */
	virtual string HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails) =0;

	size_t GetFirstIndex() const;
	size_t GetLastIndex() const;
//	bool operator<(const BlockListNode& other) const;
	/**
	 * @brief returns the number of snps contained within the block
	 */
	size_t BlockCount() const;
	int GetSelectionScore() const { return selectionScore; }
	/**
	 * @brief REturns the average minor allele frequency for the block
	 */
	double GetAverageMAF() const;

	/**
	 * @brief Return the first location on the block
	 */
	uint GetFirstLocation() const;

	/**
	 * @brief returns the last location on the bloc, 
	 */
	uint GetLastLocation() const;

	void CompleteNode(vector<Locus> &loci);
	void DeleteNode(vector<Locus> &loci);

	string GetLabel();

protected:
	BlockListNode(float criterion);
	vector<Locus> loci;						///<This is the constituent blocks only
	BlockListNode *previous;					///<The previous block in the linked list
	BlockListNode *next;						///<The next block in the linked list

	BlockList *completedList;
	BlockList *deletedList;

	uint firstIdx;								///<The start of the block
	uint lastIdx;								///<The last of the block
	float criterion;							///<This is the value we must reach before we call it a valid block member
	uint lastSnpIdxConsidered;					///<Used during splitting

	//Since we have a deleted list, we can grab one of those first before creating new
	virtual BlockListNode *NewNode(BlockListNode *previous, uint firstIdx, uint lastIdx) = 0;
	float avgMAF;								///<Average Minor Allele frequency
	size_t totalDistance;						///<Total Distance
	float blockDensity;							///<Block Density
	Locus *minAlleleFreq;						///<Minimum Minor Allele Frequency
	Locus *maxAlleleFreq;						///<Maximum Minor Allele Frequency
	int selectionScore;							///<Rank based on the fit of the validation criterion
	size_t firstLocation;
	size_t lastLocation;
	size_t previousEdge;						///<This value is used to help avoid have fake blocks that are really just pieces of other blocks

};

struct SortBlocksByLocation {
	bool operator()(const BlockListNode *left, const BlockListNode *right) {
		return left->GetFirstLocation() < right->GetFirstLocation();
	}
};

template <class T>
class BlockListHead : public BlockListNode {
public:
	BlockListHead(float criterion = 0.999);
	virtual ~BlockListHead();

	virtual void Append(uint idxA, uint idxB, vector<Locus> &loci, float value);

	void TruncateBlock(uint idxA, uint idxB);
	bool ValidateBlocks(vector<Locus> &loci, BlockFilter *filter);
	
	void SetBlockFilter(BlockFilter *filter);

	void Purge();

	BlockListNode *Clone() { assert(0); return NULL; }

	BlockList *GetValidBlocks() { return &validBlocks; }
	BlockList *GetInvalidBlocks() { return &invalidBlocks; }

	void GetBlockDensities(size_t &minBlockDensity, size_t &maxBlockDensity);

	/**
	 * @Brief starts the summary (in HTML format)
	 */
	virtual void HtmlSummaryBegin(ostream &os) {};	

	/**
	 * @brief write block details to the stream (HTML format)
	 * @param os Stream to be written to
	 * @param linkToDetails The name of the detailed plot to be imbedded in the report
	 */
	virtual void HtmlSummary(ostream &os, const char *linkToDetails) {};
		
	/**
	 * @brief Close out the block's report
	 */
	virtual void HtmlSummaryEnd(ostream &os) {};

	/**
	 * @brief Sets up the detailed hmtl report
	 * @return the name of the anchor to which the page can be linked
	 */
	virtual string HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails) { return "error"; }

	void Reset(uint first, uint last);

	BlockListNode *NewNode(BlockListNode *previous, uint firstIdx, uint lastIdx);
protected:
	BlockList validBlocks;			///List of blocks that meet criterion
	BlockList invalidBlocks;		///Lilst of blocks that fail to meet criterion
	BlockFilter *filter;
};


struct BlockSizeEval {
	bool operator()(BlockListNode* const& s1, BlockListNode* const& s2)  {
		return s2->BlockCount() < s1->BlockCount();
	}
};

struct BlNode_BlockSizeScore {
	bool operator()(const BlockListNode* a, const BlockListNode *b) const {
		if (a->BlockCount() == b->BlockCount())
			return a->GetSelectionScore() >= b->GetSelectionScore();

		return a->BlockCount() >= b->BlockCount();
	}
};


template <class T>
inline
BlockListHead<T>::BlockListHead(float criterion) : BlockListNode(criterion), filter(NULL) {
	completedList = new vector<BlockListNode *>(1024);
	deletedList = new vector<BlockListNode *>(512);
}

template <class T>
inline
void BlockListHead<T>::SetBlockFilter( BlockFilter *filter) {
	this->filter = filter;
}

template <class T>
inline
BlockListHead<T>::~BlockListHead() {
	Purge();
	if (completedList)
		delete completedList;
	if (deletedList)
		delete deletedList;
}

template <class T>
inline
BlockListNode *BlockListHead<T>::NewNode(BlockListNode *prev, uint firstIdx, uint lastIdx) {
	BlockListNode *newBlock = new T(criterion, prev, completedList, deletedList);
	newBlock->SetIndexBounds(firstIdx, lastIdx);
	return newBlock;
}

template <class T>
inline
void BlockListHead<T>::Purge() {
	BlockList::iterator cur = deletedList->begin();
	BlockList::iterator end = deletedList->end();
	
	for (; cur != end; cur++) 
		delete *cur;
	deletedList->clear();

	cur = completedList->begin();
	end = completedList->end();

	for (; cur != end; cur++) 
		delete *cur;
	completedList->clear();

	//Both of these are contained in the "completedList" vector
	validBlocks.clear();
	invalidBlocks.clear();

}

template <class T>
inline
bool BlockListHead<T>::ValidateBlocks( vector<Locus> &loci, BlockFilter *) {
	BlockList::iterator cur = completedList->begin();
	BlockList::iterator end = completedList->end();
	bool success = false;
	validBlocks.clear();
	invalidBlocks.clear();
	for (; cur != end; cur++) 
		if ((*cur)->ValidateBlocks(loci, filter)) {
			validBlocks.push_back(*cur);
			success=true;
		}
		else
			invalidBlocks.push_back(*cur);
	//sort(validBlocks.begin(), validBlocks.end(), BlNode_BlockSizeScore());
	return success;
}

template <class T>
inline
void BlockListHead<T>::GetBlockDensities(size_t &minBlockDensity, size_t &maxBlockDensity) {
	BlockList::iterator itr = completedList->begin();
	BlockList::iterator end = completedList->end();
	
	minBlockDensity = -1;
	maxBlockDensity = 0;
	while (itr != end) {
		size_t curDensity = (*itr)->GetBlockDensity();
		if (curDensity < minBlockDensity) 
			minBlockDensity = curDensity;
		if (curDensity > maxBlockDensity)
			maxBlockDensity = curDensity;
		itr++;
	}
}

template <class T>
inline
void BlockListHead<T>::Reset(uint first, uint last) {
	//cout<<"Resetting head: "<<firstIdx<<" - "<<lastIdx<<" : "<<first<<" - "<<last<<"\n";
	firstIdx = first;
	lastIdx = last;
	Purge();
	completedList->clear();
	deletedList->clear();
	validBlocks.clear();
	invalidBlocks.clear();	
}

template <class T>
inline
void BlockListHead<T>::Append(uint idxA, uint idxB, vector<Locus> &loci, float value) {
	if (next == NULL) {
		next = NewNode(this, firstIdx, lastIdx);
	} 
	next->Append(idxA, idxB, loci, value);
}

template <class T>
inline
void BlockListHead<T>::TruncateBlock(uint idxA, uint idxB) {
	if (next == NULL) 
		next = NewNode(this, firstIdx, lastIdx);
	next->TruncateBlock(idxA, idxB);
}


}

#endif
