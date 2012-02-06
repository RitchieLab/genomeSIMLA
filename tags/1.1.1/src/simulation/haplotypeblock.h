//
// C++ Interface: haplotypeblock
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONHAPLOTYPEBLOCK_H
#define SIMULATIONHAPLOTYPEBLOCK_H
#include "locus.h"
#include <iostream>
#include <map>
#include "simulation.h"

namespace Simulation {

using namespace std;

/**
	@brief Represents details associated with a haplotype block
	@note Create a block with first and last snp- and grow it by appending loci
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class HaplotypeBlock{
public:

	/**
	 * @brief Help keep track of the loci within the original array
	 */
	struct LocusWrapper {
		public:
		LocusWrapper(size_t idx, Locus *loc) {
			locus = loc;
			index = idx;
		}

		Locus *locus;
		size_t index;
	};

	/**
	 * @brief functor to help sort Haplotype blocks based on MAF
	 */
	struct LTMAF {
		bool operator()(const LocusWrapper& a, const LocusWrapper& b) const {
			return a.locus->GetMinAlleleFreq() < b.locus->GetMinAlleleFreq();
		}
	};

	/**
	 * @brief Functor to help sort haplotype blocks based on position
	 */
	struct LTPos {
		bool operator()(const LocusWrapper& a, const LocusWrapper& b) const {
			return a.locus->GetLocation() < b.locus->GetLocation();
		}
	};

	/**
	 * @brief Initialize a block with two loci and their respective indices
	 */
	HaplotypeBlock(size_t idx1, Locus *loc1, size_t idx2, Locus* loc2);

	~HaplotypeBlock();

	/**
	 * @brief Append a new locus to the block
	 */
	void Append(size_t idx1, Locus *loc1);
	
	/**
	 * @brief add a probability structure between two snps (for reporting purposes)
	 */
	void AppendProbabilities(const char *snp1, const char *snp2, Probabilities &prob1);

	/**
	 * @brief generate summary report
	 */
	void Report(ostream& os);

	/**
	 * @brief Generate a more detailed report of the block contents
	 */
	void DetailedReport(ostream &os);

	/**
	 * @brief returns the number of snps contained within the block
	 */
	size_t BlockCount() const;

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

	size_t GetFirstIndex() const;
	size_t GetLastIndex() const;

	/**
	 * @brief Returns the block density
	 */
	size_t GetBlockDensity() const;
	
	/**
	 * @brief This is used for sorting in STD sort functions
	 * @note We want it to work out so the best block is the first in the vector
	 */
	bool operator<(const HaplotypeBlock& other) const;
	uint firstLocusIndex;					
	uint lastLocusIndex;					
	string firstSnpLabel;
	string lastSnpLabel;

	/**
	 * @Brief starts the summary (in HTML format)
	 */
	void HtmlSummaryBegin(ostream &os);	

	/**
	 * @brief write block details to the stream (HTML format)
	 * @param os Stream to be written to
	 * @param linkToDetails The name of the detailed plot to be imbedded in the report
	 */
	void HtmlSummary(ostream &os, const char *linkToDetails);
		
	/**
	 * @brief Close out the block's report
	 */
	void HtmlSummaryEnd(ostream &os);

	/**
	 * @brief Sets up the detailed hmtl report
	 * @return the name of the anchor to which the page can be linked
	 */
	string HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails);
	static uint maxImageWidth;				///<Used to size the images on the report pages
protected:
	vector<LocusWrapper> loci;				///<Vector of loci
	Locus *minAlleleFreq;					///<This is the locus with the lowest MAF (non-zero)
	Locus *maxAlleleFreq;					///<This is the locus with the highest MAF
	double totalMAF;						///<This is the sum of minor allele frequencies (used for average)

	uint lastIndex;

	/**
	 * @brief Keep track of the various probabilities associated with each piece of the block
	 */
	map<string, Probabilities> probLookup;


//For HTML columns
	static uint rowsWritten;

};

struct BlockSizeEval {
	bool operator()(HaplotypeBlock*& s1, HaplotypeBlock*& s2) {
		return s2->BlockCount() < s1->BlockCount();
	}
};


}

#endif
