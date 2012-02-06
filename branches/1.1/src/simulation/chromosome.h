//Chromosome.h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// 
//
// Stores allele information for a chromosome.  
// Static functions and data allow for common genome to be
// represented across the entire set of chromosomes.
//
/////////////////////////////////////////////////////////////////////


#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include "boost/dynamic_bitset.hpp"

#include "utility/random.h"
#include "locus.h"
#include "locusassociation.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "utility/types.h"
#include "utility/rbtree.h"
#include <map>
#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif

struct uintLT {
	int operator()(const uint l, const uint r) const {
		if (l<r) return -1;
		if (l>r) return 1;
		return 0;
	}
};

struct floatLT { 
	int operator()(const float l, const float r) const {
		if (l<r) return -1;
		if (l>r) return 1;
		return 0;
	}
};
typedef Utility::RBTree<float, int, floatLT> RecDistType;
typedef Utility::RBTreeNode<float, int, floatLT> RecDistNodeType;

typedef Utility::RBTree<size_t, std::string> AllSrcType;
typedef Utility::RBTreeNode<size_t, std::string> AllSrcNodeType;

namespace Simulation {

class MinChrom;



typedef std::vector<Locus> LocusArray;

/**
 * @brief Manage the chromsomal information. 
 * Each chromosome depends on a predefined array of loci. Two chromosomes can be crossed to simulate
 * sexual reproduction
 */
class Chromosome{
friend class ChromosomeTest;
  public:
	/**
	 * @brief Construction
	 * @param loci The vector of loci objects
	 * @param pLamdba The poisson lambda. 
	 * @param recombIndexLookup Reverse lookup of indices by location
	 * @note We don't want to calculate lambda ourselves, since this assignment will happen over and over
	 **/
	Chromosome(LocusArray *loci, double pLambda, RecDistType *recombIndexLookup);
	/**
	 * @brief Constructor...should only be used by STL
	 */
	Chromosome();

	/**
	 * @brief Copy Constructor. This is a complete copy...down to the xo-events
	 */
	Chromosome(const Chromosome& other);

	void Distort(const double& freq);

	/**
	 * @brief Special constructor to allow delayed cross over to occur
	 */
	Chromosome(const Chromosome &mother, const Chromosome &father);
	~Chromosome() {
		if (xoEvents)	delete xoEvents;
	}
	void ReferenceDistance(std::vector<double>& distances, int diseaseLocus);
	void XOCount(std::vector<int>& counts, int diseaseLocus);
	void EvaluateRecombinants(std::vector<int>& r, int diseaseLocus);
	/**	
 	 * @brief Returns a reference to the bit at position (locusIndex)
	 * @param locusIndex The position within the chromosome
	 */
    boost::dynamic_bitset<>::reference operator [] (uint locusIndex) {
		assert(locusIndex < chrom.size());
	
		if (hasCompleteGenotypes)
			return chrom[locusIndex];
		else {
//assert(0);
			//We need to figure out which phase we the genotype is at and ask our parent chromosome
			if (GetPhase(locusIndex))
				return (*phasedChrom)[locusIndex];
			else
				return chrom[locusIndex];
		}
	}		

	/**
	 * @brief Returns the phase of a given locus when using delayed crossing
	 * @return 0/1 for whether it's the other chromosome (phase 0 has already been copied)
	 */
	bool GetPhase(int locusIndex) {
		//assert(!hasCompleteGenotypes);
		bool thePhase = false;
		Utility::RBTreeNode<size_t, bool> *node = genotypePhase.FindNearestMin(locusIndex);
		if (node) {
			thePhase = node->GetData();
	//		std::cout<<"GetPhase("<<locusIndex<<") = Idx: "<<node->GetKey()<<" : "<<node->GetData()<<"\n";
		}
		return thePhase;
	}

	void RenderPhase(std::ostream& os, uint count) {
		if (count > chrom.size()) 
			count = chrom.size();
		for (uint i=0; i<count; i++) 
			os<<GetPhase(i)<<" ";
	}
	
	Chromosome &operator=(const Chromosome& other);
	Chromosome &operator=(MinChrom& other);
	void ResolveTo(const Chromosome &other);
	/**
	 * @brief Generate a new set of cross over events. 
	 * @return Returns a pointer to an array of indices. 
	 * @note Each value in the array indicates a locus whos phase is different from the one before.
	 */
	std::vector<size_t> *BuildXOEvents() const;
	/**
	 * @brief Represent sexual reproduction
	 * @param secondChrom The other "parent" involved with the crossing
	 * @return New chromosome. This should be deleted when it is time
	 */
    //Chromosome Cross(Chromosome &secondChrom);
    Chromosome Cross(Chromosome &secondChrom, size_t &crossoverEvents);

	void ResetToMinimal(std::vector<uint> &modelLoci);
	void ResolveGenotypes();

	/**
	 * @brief streaming operator
	 */
    friend std::ostream & operator << (std::ostream & os, Chromosome & chrom);
  
	/**
	 * @brief Write the chromosome in ped format
	 */
	void WritePedFormat(std::ostream& os, float threshold, uint first, uint last);

	void WriteXOPoints(std::ostream& os);
	
	/**
	 * @brief Returns number of loci assicated with the local "pool"
	 */
    size_t LociCount() const {if (loci) return loci->size(); else return 0;}

	Chromosome CrossImmediate(Chromosome &secondChrom, size_t &crossoverEvents);
	Chromosome DLCross(Chromosome &secondChrom, size_t &crossoverEvents);
	/**
	 * @brief return a reference to the local locus at the desired position
	 * @param index position in the locus array
	 */
    Locus& GetLocus(uint index){return (*loci)[index];}

	/**
	 * @brief returns the allelic value associated with index, locIdx
	 */
	uint At(int locIdx) { assert((uint)locIdx < chrom.size()); return (*this)[locIdx]; }

	/**
	 * @brief Eventually, we'll want to assign specific random number generators to each part of the system
	 */
	//static Utility::Random *generator;

	/**
	 * @brief Intialize a given locus
	 */
	void InitLocus(int locusIndex);

	/**
	 * @brief Initialize all loci 	
	 */
	void InitLoci(LocusAssociationGrid *grid);
	void Invert();
	void InitLoci(std::vector<int>& alleles);

	void ExpandModelLoci(std::vector<uint> &modelLoci);

	void WriteBinary(std::ostream& os);
	void ReadBinary(std::istream& is, size_t lociCount, std::vector<uint> *modelLoci, bool retainChrom);
//	void ReadBinary( std::ifstream *file, uint length);
	
	void ShowGenotypes(uint count, char div);

	/**
	 * @brief Draw an event from the poisson distribution based on the chromosome's total rec. frequency
	 * @param lambda This is the mean number of recombinations to occur over the length of the chromosome
	 * @return number of recombination events that might occur over the entire chromosome
	 */
	size_t PoissonEventCount() const;

	/**
	 * @brief force the value upon the data
	 */
	void SetValue(uint locus, uint value);
	void InitAlleleSource(const char *source_id);

	/**
	 * @brief Returns the number of SNPs that have a common source between the two chromosomes
	 */
	int EvaluateKinship(Chromosome& other);
  protected:
	Utility::Random *generator;
    boost::dynamic_bitset<> chrom;				///<The allelic data
	LocusArray *loci;							///<Describes various loci
	double poissonLambda;						///<Lambda to be used by poisson distribution
	RecDistType *recombIndexLookup;				///<Quickly lookup the index at a given location
	std::vector<size_t> *xoEvents;				///<Index of cross over events
	bool hasCompleteGenotypes;					///<Indicates whether the chromosomal data has been completely flushed out
	Utility::RBTree<size_t, bool> genotypePhase;
	Utility::RBTree<size_t, std::string> alleleSource;
	Chromosome *phasedChrom;
	Chromosome *unphasedChrom;					///<Only used when we cross over from an empty pool (which needs to be loaded)
	bool phase;									///<Indicates which phase we start on (true means copy from phasedchrom)
	//boost::dynamic_bitset<> modelOnlyChrom;		///<Keep up with only the loci associated with a given model
	std::map<uint, bool> modelOnlyChrom;

	void ExtractAlleleSource();
};




#ifdef CPPUNIT
class ChromosomeTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( ChromosomeTest );
//	CPPUNIT_TEST( TestKinship );
	CPPUNIT_TEST_SUITE_END();
public:
    ChromosomeTest();
    ~ChromosomeTest();


	void setUp();
	void tearDown();
	
	void TestKinship();
protected:
	LocusArray loci;
	Chromosome *ch1;		//Founder
	Chromosome *ch2;		//Founder

	Chromosome *ch3;		//Sib
	Chromosome *ex1;		//External 

	Chromosome *ch4;		//Sib
	Chromosome *ex2;		//External 

	Chromosome *ch5;		//Grandchild
	Chromosome *ex3;		//External 

	Chromosome *ch6;		//Grandchild
	Chromosome *ex4;		//External 

	Chromosome *ch7;		//Great Grandchild
	Chromosome *ch8;		//Great Grandchild
};
#endif

}


#endif
