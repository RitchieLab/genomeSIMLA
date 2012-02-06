//
// C++ Interface: ldcalculator
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDCALCULATOR_H
#define SIMULATION_VISUALIZATIONLDCALCULATOR_H

#include "blocklistnodefourgammetes.h"
#include "allelesource.h"
#include "utility/array2d.h"

#define dotA  P->Freq1()
#define dota  P->Freq2()
#define dotB  Q->Freq1()
#define dotb  Q->Freq2()

namespace Simulation {

namespace Visualization {


/**
 * @brief Storage for a single set of LD results between two SNPs
 */
struct LdResult {
public:
	LdResult();

	/**
	 * @Brief initialize with frequency counts and pointers to the relevant loci
	 */
	LdResult(Locus*p, Locus*q, int cab, int caB, int cAb, int cAB);
	Locus *P;
	Locus *Q;

	float DPrime(float pab, float paB, float pAb, float pAB);
	float RSquared(float pab, float paB, float pAb, float pAB);
	float LOD(int cab, int caB, int cAb, int cAB, float pab, float paB, float pAb, float pAB);
	bool CheckMarginals(Locus *p, Locus* q, float pab, float paB, float pAb, float pAB);
	float dprime;
	float rsquared;
	float lod;	
};

extern unsigned int maxSnpDistance;			///<Cutoff for ld calculations


/**
 * 	@brief Calculates all LD associate with the first locus and all subsequent loci
 * 
 *	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>

	LdCalculator will calculate ALL Ld values for the SNPs in ld with firstLocus up to (not including) lastLocus;
	This is a single diagonal on the inverted LD Pyramid.

	 * Usage			- The pool instantiates one or more copies of LdCalculator, providing the desired 
	 *						chromosomes. It is possible to instantiate multiple copies for the same 
	 * 						region and feed each one different chromosomes and then merge them all prior 
	 * 						to LD extraction. The Haplotype block object is used during ld calculation to
	 * 						evaluate the different blocks encountered. The results represent a single 
	 * 						diagonal and only contain nodes for which LD values were calculated. 
	
	 * Initialization	- The structure is initialized by the region (in vector iterators) within a
	 * 						vector of Locus pointers. The calculations are only done during CalculateFrequencies
	 * 						though. 

	 * Reuse			- The structure can be reused for different regions (i.e. dropping the first locus
	 * 						from the calculation. Simply recall *Init* and *CalculateFrequencies* and += (as necessary)
	
	 * Threading - 		- LdCalculator Objects can be summed *AFTER* CalculateFrequencies have been called, as long
	 * 						as they represent the same Locus region. If this is the case, a+=b will result in a 
	 * 						being valid for LD associated with individuals seen by both a and b. Once all have
	 * 						been "added" together, GetLdResults() will be appropriate for all added instances
	
	 * Memory			- Memory consumption is non-trivial. Once initialized, an array of floats sized 
	 * 						4 times the size of the number of loci is instantiated. These values are 
	 * 						deleted prior to each subsequent run and do not accumulate. During the final step
	 * 						N LdResults will be appended to the vector<LdResults>, during which time the structure
	 * 						continues to grow. 
	 * 
     * 					- In addition to general memory requirements, the size of the result vector will grow
	 * 						with each iteration. To maximize efficiency, one should reserve() a reasonable number
	 * 						in the initial vector in order to avoid the vector's recopying, which could slow things
	 * 						down. one solution to this problem would be to use the formulat: reserve(N*N/2 + n/2)
	 * 						prior to the first call. This would size the number of cells equal to the largest
	 * 						possible number necessary and avoid any reallocation.
 
 */
template <class T>
class LdCalculator {
public:
	LdCalculator(typename vector<AlleleSource<T> *>::iterator firstChrom, typename vector<AlleleSource<T>* >::iterator endChrom, vector<Locus> *loci, BlockListHead<BlockListNodeFourGammetes> *haplotypeHead);

	/**
	 * @Brief Establish the boundaries between the first and last snp in the plot 
	 * @param start the first locus to be considered
	 * @param end the end clause...LD with this snp will not be evaluated (because it might be .end()
	 */
	bool Init(int start, int end);

	/**
	 * @Brief Perform the actual frequency calculation based on the individuals found in the iterators
	 */
	void CalculateFrequencies();

	/**
	 * @Brief sum together the frequencies, so that you can put multiple threads together
	 */
	LdCalculator& operator+=(const LdCalculator& other);
	
	/**
	 * @brief Obtain the actual LD statistics based on the frequency data found in the local object
	 * @return Returns the list of haplotype blocks that were encountered
	 */
	void GetLDResults(vector<LdResult>& results);

protected:
	///These are the chromosomes used to seed the frequency data
	typename vector<AlleleSource<T> *>::iterator firstChrom;
	typename vector<AlleleSource<T> *>::iterator lastChrom;

	///These are the loci of interest. The snp at firstLocus is used to compare with all subsequent ones
	int firstLocus;
	int lastLocus;

	int locusCount;			///<The number of loci to be processed
	Array2D<int> freq;				///<Frequencies of each pairing

	/**	
	 * @brief aide in the production of haplotype blocks
	 */
	BlockListHead<BlockListNodeFourGammetes> *haplotypeHead;

	vector<Locus> *loci;
};


template <class T>
inline
LdCalculator<T>::LdCalculator(typename vector<AlleleSource<T> *>::iterator firstChrom, 
			typename vector<AlleleSource<T> *>::iterator lastChrom, vector<Locus>* loci, BlockListHead<BlockListNodeFourGammetes> *haplotypeHead) : firstChrom(firstChrom), lastChrom(lastChrom), firstLocus(0), lastLocus(0), locusCount(0), loci(loci), haplotypeHead(haplotypeHead) { }

template <class T>
inline
bool LdCalculator<T>::Init(int first, int last) {
	firstLocus = first;
	lastLocus = first+1;
	int idx  = first;
	locusCount=0;
	
	//An empty vector?
	if (first == last) 
		return false;
	
	Locus *loc = &loci->at(first);
	uint maxLoc = loc->GetLocation() + maxSnpDistance;
	bool maxFound=false;

	for (; lastLocus<=last && !maxFound; ) {
		loc = &loci->at(lastLocus);
		maxFound = maxLoc < loc->GetLocation();
		if (maxFound)
			lastLocus--;
		else
			lastLocus++;
	}
	if (lastLocus > last)
		lastLocus = last;
	locusCount = lastLocus - firstLocus;

	freq = Array2D<int>(locusCount, 4);
	//We don't evaluate the last node, if there are none between the first and last, then it's an empty set
	if (locusCount == 0)
		return false;
	return true;
}

template <class T>
inline
//Simple, no parameter call...perfect for being called from a thread.....just call one LdCalculator per thread
void LdCalculator<T>::CalculateFrequencies() {
	Locus &A = loci->at(firstLocus);
	int idxA = A.GetID();

	//Iterator over each individual
	typename vector<AlleleSource<T> *>::iterator chrom = firstChrom;
	while (chrom != lastChrom) {
		AlleleSource<T> *c = *chrom;
		int aVal = c->At(idxA) *2;
		
		int idx = 0;
		int lastValue = 1;
		for (int i=firstLocus+1; lastValue > 0 &&  i<=lastLocus; i++) {
			//Attempt to short circuit the calculations
			freq(idx++, aVal+c->At(loci->at(i).GetID()))++;
		}
		chrom++;
	}
}

template <class T>
inline
void LdCalculator<T>::GetLDResults(vector<LdResult>& results) {
	//Try to avoid reallocating the vector more than once...
	int curSize = results.size();
	//This should be OK, even if you've already reserved enough to hold the entire structure...
	results.reserve(locusCount + curSize);

	Locus &A = loci->at(firstLocus);
	int idxA = A.GetID();
	int idx = 0;
	for (int i=firstLocus+1; i<=lastLocus; i++) {
		Locus &B = loci->at(i);
		LdResult result(&A, &B, freq(idx,0), freq(idx,1), freq(idx,2), freq(idx,3));
//assert(result.dprime > 0.0);
		results.push_back(result);
		haplotypeHead->Append(firstLocus, i, *loci, result.dprime);
		idx++;
	}
	if (idx == 0) 
		haplotypeHead->Append(firstLocus, firstLocus+1, *loci, 0);
	haplotypeHead->TruncateBlock( firstLocus, lastLocus+1);
}

template <class T>
inline
LdCalculator<T>& LdCalculator<T>::operator+=(const LdCalculator<T>& other) {
	assert(locusCount == other.locusCount);
	int count = locusCount *4;
	freq=other.freq;
	return *this;
}



}

}

#endif
