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
#include <iostream>
#include <math.h>

namespace Simulation {

class Chromosome;

typedef std::vector<Locus> LocusArray;

/**
 * @brief Manage the chromsomal information. 
 * Each chromosome depends on a predefined array of loci. Two chromosomes can be crossed to simulate
 * sexual reproduction
 */
class Chromosome{
  public:
	/**
	 * @brief Construction
	 * @param loci The vector of loci objects
	 **/
    Chromosome(LocusArray *loci) : loci(loci)	{ 	
		generator = &Utility::Random::globalGenerator;
		chrom.resize(loci->size()); 
	}

	Chromosome() : loci(NULL) {
		generator = &Utility::Random::globalGenerator;
	}

	/**	
 	 * @brief Returns a reference to the bit at position (locusIndex)
	 * @param locusIndex The position within the chromosome
	 */
    boost::dynamic_bitset<>::reference operator [] (int locusIndex){return chrom[locusIndex];}
    
	/**
	 * @brief Represent sexual reproduction
	 * @param secondChrom The other "parent" involved with the crossing
	 * @return New chromosome. This should be deleted when it is time
	 */
    Chromosome Cross(Chromosome &secondChrom);
    
	/**
	 * @brief streaming operator
	 */
    friend std::ostream & operator << (std::ostream & os, Chromosome & chrom);
  
	/**
	 * @brief Write the chromosome in ped format
	 */
	void WritePedFormat(std::ostream& os, uint first, uint last);
	
	/**
	 * @brief Returns number of loci assicated with the local "pool"
	 */
    unsigned int LociCount(){return loci->size();}

	/**
	 * @brief return a reference to the local locus at the desired position
	 * @param index position in the locus array
	 */
    Locus& GetLocus(uint index){return (*loci)[index];}

	/**
	 * @brief returns the allelic value associated with index, locIdx
	 */
	uint At(int locIdx) { return chrom[locIdx]; }

	/**
	 * @brief Eventually, we'll want to assign specific random number generators to each part of the system
	 */
	static Utility::Random *generator;

	/**
	 * @brief Intialize a given locus
	 */
	void InitLocus(int locusIndex);

	/**
	 * @brief Initialize all loci 	
	 */
	void InitLoci();

	/**
	 * @brief force the value upon the data
	 */
	void SetValue(uint locus, uint value);
  private:
    boost::dynamic_bitset<> chrom;				///<The allelic data
	LocusArray *loci;							///<Describes various loci
};








}


#endif
