//
// C++ Interface: locusselection
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLOCUSSELECTION_H
#define SIMULATIONLOCUSSELECTION_H

#include <vector>
#include <sstream>
#include "utility/types.h"
#include "locus.h"
#include "blocklistnode.h"
#include "utility/rbtree.h"

namespace Simulation {

namespace Visualization {

using namespace std;


/**
 * @brief Helper class that determines if a value fits within the range and returns a value from 0.0 to 1.0 indicating how far away from the goal it is
 */
template <class T>
struct LinearRange {
	LinearRange(const LinearRange& other) : target(other.target), min(other.min), 
			max(other.max), maxBound(other.maxBound), minBound(other.minBound) { }
	LinearRange(T target, T min, T max) : target(target), min(min), max(max) {
		
		maxBound = max * (1+1e-32);
		minBound = min * (1-1e-32);
	}
	
	float Evaluate(T value);
	
	T target;
	T min;
	T max;
 	double maxBound;
	double minBound;

	string GetDescription() {
		stringstream ss;
		ss<<target<<" : "<<min<<" - "<<max;
		return ss.str();
	}
};

/**
 * @brief Returns true or false if the value falls within the range
 */
template <class T>
struct SimpleRange {
public:
	SimpleRange(const SimpleRange& other) : min(other.min), max(other.max) { }
	SimpleRange(T min, T max) : min(min), max(max) { }
	bool Evaluate(T value) { return (value >= min && value <= max);  }

	T min; 	
	T max;

	string GetDescription() {
		stringstream ss;
		ss<<min<<" - "<<max;
		return ss.str();
	}

};

struct LocusReportEntry {
	float score;						///<The score for this entry
	Locus locus;						///<The locus of interest
	string blockDescriptionFile;		///<File where the block is described
	string chromosomeID;				///<For reporting purposes
	~LocusReportEntry();				
	
	/**
	 * @brief Compares the scores ( > )
	 */
	bool operator()(const LocusReportEntry& l, const LocusReportEntry& r) { return l.score < r.score; }

	bool operator<(const LocusReportEntry& other) { return score < other.score; }
	BlockList associatedBlocks;
		
	/**
	 * @brief Writes html report to the stream
	 * @param report The stream to be written to
	 * @param useAltColor Allow for alternate row shading
	 */
	void WriteHtmlReport(ostream& report, bool useAltColor);

	/**
	 * @brief Write the header for the report to the stream
	 * @param file Stream to be written to
	 */
	void WriteHtmlHeader(ostream &file);

	void WriteToText(ostream &file);
};

/**
	@brief Used to rank loci for final reporting (as well as blocks)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LocusSelection{
public:
	/**
	 * @brief required by map. Please use the other one
	 */
	LocusSelection()  : mafRange(0.2, 0.15, 0.25), blockSize(3, 2, 5) { }

	LocusSelection(const LocusSelection& other) : mafRange(other.mafRange), blockSize(other.blockSize), ranges(other.ranges), contents(other.contents), label(other.label), description(other.description) {  }
    LocusSelection(float maf, float min, float max, uint blck, uint bMin, uint bMax) : mafRange(maf, min, max), blockSize(blck, bMin, bMax) {}
    ~LocusSelection() {	}

	/**
	 * @brief Evaluate a given locus for the score based on local criterion
	 * @param loc Locus of interest
	 * @param blocks vector of blocks where loc is present
	 */
	float Evaluate(LocusReportEntry* locReport);
	void WriteToText(ostream& report);
	/**
	 * @brief Copy all range information from the other selection object
	 */
	void ExtractRanges(LocusSelection &other);

	/**
	 * @brief Evaluates a locus for a given range
	 */
	bool EvaluateRange(Locus &loc);

	/**
	 * @brief Just a quick check to make sure the locus falls within the range
	 */
	//bool WithinRange(Locus &loc);

	/**
	 * @brief Add two loci as starting/stopping points. If any points exist between the two points, they will be removed (or expanded)
	 */
	bool AddRange(Locus& start, Locus& stop);

	bool SetRange(int idx, Locus& start, Locus &stop);


	static float mafWeight;					///<weigth used to evaluate deviation from target MAF
	static float blockSizeWeight;			///<Weight used to evaulate deviation from target block size
	static uint maxReportedEntries;			///<controls how many Loci are reported

	/**
	 * @brief Write the report to the stream
	 */
	void WriteHtmlReport(ostream& report);

	/**
	 * @brief Sets the description which can be a multiword phrase
	 */
	void SetDescription(const char *desc);

	/**
	 * @brief Set the label for the selector object.
	 * The label is used for quick identification
	 */
	void SetLabel(const char *label);

	/**
	 * @brief Returns the label
	 */
	string GetLabel() { return label; }

	/**
	 * @brief Returns the description	
	 */
	string GetDescription() { return description;  }

	/**
	 * @Brief Resets the contents associated with the selector.
	 * Use in between drops
	 */
	void Reset();

	string GetDescMAF();
	string GetDescBlockSize();
	string GetDescRanges();

	void SetBlockSizeRange(LinearRange<uint> range);
	void SetMAFRange(LinearRange<float> range);

	LinearRange<float> &GetMafRange() 	{ return mafRange; }
	LinearRange<uint> &GetBlockSizeRange() { return blockSize; }

	SimpleRange<Locus> &GetRange(size_t idx) { assert(idx < ranges.size()); return ranges[idx]; }
	size_t GetRangeCount() { return ranges.size(); }
	bool RemoveRange(size_t idx);

	string GetConfiguration(const char *delim);
protected:
	LinearRange<float> mafRange;				///<Range of tolerable MAFs
	LinearRange<uint> blockSize;				///<Range of tolerable Block Sizes
	std::vector<SimpleRange<Locus> > ranges;	///<All the ranges that can be considered
	
	/**
	 * @brief The entries that match our needs, sorted by their score
	 */
	Utility::RBTree<float, LocusReportEntry* > contents;

	string label;								///<The label for quick identification
	string description;							///<Phrase that describes this selector object
	
};


template <class T>
inline
float LinearRange<T>::Evaluate(T value) {
	float score = 0.0;

	if (value == target)
		return 1.0;
	else if (value < min || value > max)
		return 0.0;

	//We don't have to worry about the denominator here, due to the first 
	//evaluation (value == target)
	if (value < target) 
		score = (double)(value - minBound + 1)/(double)(target - minBound + 1);
	else 
		score = (double)(maxBound - value + 1) / (double)(maxBound - target + 1);

	return score;
}

}

}

#endif
