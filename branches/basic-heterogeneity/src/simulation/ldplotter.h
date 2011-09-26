//
// C++ Interface: ldplotter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLDPLOTTER_H
#define SIMULATIONLDPLOTTER_H
//#include <pngwriter.h>
#include <vector>

#include "locus.h"
#include "chromosome.h"
#include "ldwriter.h"
//#include "haplotypeblock.h"
#include "blocklistnodefourgammetes.h"
#include "simulation.h"
namespace Simulation {

namespace Visualization {

using namespace std;

/**
 * @brief Cache entry for single ld statistic
 */
struct LDStatistics {
	Locus *A;					///<First Locus
	Locus *B;					///<Second Locus
	float dprime;				///<Dprime value
	float lod;					///<LOD 
	float rsquared;				///<R-Squared Value
	bool deadSpace;				///<Used to indicate we have a blank spot
	LDStatistics(bool isDead = true) : A(NULL), B(NULL), dprime(0.0), lod(0.0), rsquared(0.0), deadSpace(isDead) { }
	LDStatistics(Locus *A, Locus *B, float dprime, float lod, float rsq) : A(A), B(B), dprime(dprime), lod(lod), rsquared(rsq), deadSpace(false) { }
};
/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
	@brief Used to help render the various LD Plots
	@Note This class is responsible for actually calculating the LD values and sending them to the various writer objects
*/
class LdPlotter{
public:

	/**
	 * @brief Construction
	 * @param thresh Trash- this needs to be removed
	 * @param pool The pool of chromosomes required to calculating frequency counts
	 * @param loci The array of loci. These should already have their frequencies calculated. 
	 */
    LdPlotter(float thresh, vector<Chromosome>& pool, vector<Locus>& loci);

    ~LdPlotter();

	/**
	 * @brief Initialize the plotter and return the vector of Haplotype blocks
	 */
	//vector<HaplotypeBlock *>* Init();
	BlockList *Init();
	
	/**
	 * @brief Determine the maximum depth of a given range of SNPs will require to plot them onto the various LD plots
	 */
	void CountRC(size_t first, size_t last, size_t &rows, size_t &cols);
	void CountRC(size_t &rows, size_t &cols);

	/**
	 * @brief Produce a DPrime ld plot
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @param imgHeight Returns the actual height of the chart
	 * @param imgWidth Returns the actual width of the chart
	 * @param ldBufferSize The number of snps on either side of the block that will show up in the plot
	 * @param Lbl Text label used as part of the final filename
	 * @return Actual filename produced
	 */
	string WriteLdDPrime(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl);
	//string WriteLdDPrime(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl);

	/**
	 * @brief Produce a ld plot based on rsquared
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @param imgHeight Returns the actual height of the chart
	 * @param imgWidth Returns the actual width of the chart
	 * @param ldBufferSize The number of snps on either side of the block that will show up in the plot
	 * @param Lbl Text label used as part of the final filename
	 * @return Actual filename produced
	 */
	string WriteLdRSquared(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl);
	//string WriteLdRSquared(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl);

	/**
	 * @brief Write the LD values to a plain text file	 
	 * @param filename This is the first portion of the filename. The actual filename will be returned
	 * @param block The block around which the chart will be centered
	 * @return Actual filename produced
	 */
	string WriteLdReport(const char *filename, BlockListNode *block);
	//string WriteLdReport(const char *filename, HaplotypeBlock *block);

	/**
	 * @brief Generic LD production routine- ld values for region surrounding the block are calculated (or reused)
	 * @param filename Basis for the actual filename to be used
	 * @param block The block around which the plot will be built
	 * @param writer The object that will actually render the desired details (text or graphical)
	 * @param maxSnpDepth The most snps a single snp can be paired with
	 * @return Actual filename used
	 */
	string WriteLdData(const char *filename, BlockListNode *block, LdWriter* writer, uint maxSnpDepth);
	//string WriteLdData(const char *filename, HaplotypeBlock *block, LdWriter* writer, uint maxSnpDepth);

	/**
	 * @brief Generic LD production routine- ld values for region surrounding the block are calculated (or reused)
	 * @param filename Basis for the actual filename to be used
	 * @param first The index of the first snp to be rendered
	 * @param last The index of the last snp to be rendered
	 * @param writer The object that will actually render the desired details (text or graphical)
	 * @param maxSnpDepth The most snps a single snp can be paired with
	 * @return Actual filename used
	 */
	string WriteLdData(const char *filename, uint first, uint last, LdWriter* writer, uint maxSnpDepth);

	static uint maxSnpDistance;			///<The max number of basepairs allowed for LD Consideration

	
	/**
	 * @brief Set a label for the chart
	 */
	void SetLabel(const char *lbl);

	void CalculateLD();
	
	/**
	 * @brief calculates LD using poisson distribution. 
	 * @param filename The LD values will be written (for drawing LD/distance plots) 
	 */
	void CalculateLD2(const char *filename);

protected:		
	/**
	 * @brief Calculate Allele frequencies for the given loci
	 */
	void CalculateAlleleFrequencies(vector<Locus>& loci);

	/**
	 * @brief Count Haplotype frequencies
	 * @param A Index of the first snp
	 * @param B Index of second snp
	 * @param p The probabilities assocaited with the two snps
	 */
	void CountFrequencies(uint A, uint B, Probabilities &p);
	void CountFrequencies(uint A, uint B1, uint B2, Probabilities* &p);
	void CountFrequencies(uint A, uint B1, uint B2, Probabilities* &p, uint firstIndividual, uint totalInd);
	Probabilities *ThreadedFrequencyCount(uint A, uint B1, uint B2, vector<pthread_t>& threads);

	static void *CountFrequencies(void *);
	/**
	 * @brief Determines the haplotype blocks for the loci within the region
	 */
	//uint BuildHaplotypeBlocks();

	/**
	 * @brief Vector of haplotype blocks that will be used by client objects to perform certain tasks
	 */
	//vector<HaplotypeBlock *> haplotypes;
	BlockList *haplotypes;
	

	vector<Chromosome> pool;			///<Required for calculating haplotype frequencies
	vector<Locus> loci;					///<The loci themselves
	float threshold;					///<Minimum MAF

	size_t minBlockDensity;				///<Min. Snp density seen for the region, used for calculating the SNP density bars
	size_t maxBlockDensity;				///<Max. Snp density seen for the region, used for calculating the SNP density bars
	string label;						///<Label that goes onto the plot
	LDStatistics *ldStatistics;			///<Array of statistics
	size_t *ldStart;					///<The index into the array for a given snp (this array is 1 entry 
										///<larger than the number of legal sub arrays

	BlockListHead<BlockListNodeFourGammetes> haplotypeHead;
	

};
}
}

#endif
