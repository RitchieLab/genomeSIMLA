// Chromosome.cpp

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

#include "chromosome.h"

namespace Simulation {

Utility::Random *Chromosome::generator = NULL;


// Use: Initializes specified locus according to allele
//      frequencies
// Arg: locusIndex -- index of locus to initialize
// ret: none
void Chromosome::InitLocus(int locusIndex){
	assert(generator);
	if(generator->drand()<=(*loci)[locusIndex].GetAlleleFreq(1))
    	chrom[locusIndex]=0;
	else
    	chrom[locusIndex]=1;
}


// Use: Initializes all loci in chromosome according to allele
//      frequencies
// Arg: none
// ret: none
void Chromosome::InitLoci(){
	unsigned int totalLoci = Chromosome::loci->size();
	for(unsigned int currLocus=0; currLocus < totalLoci; currLocus++){
 		if(generator->drand()<=(*loci)[currLocus].GetAlleleFreq(1))
     		chrom[currLocus]=0;
    	else
      		chrom[currLocus]=1;
  	}
}


// Use: Produces new chromosome by crossing this chromosome with
//      a second one.  
// Arg: secondChrom - chromosome to cross
// Ret: Pointer to new chromosome
Chromosome Chromosome::Cross(Chromosome &secondChrom){
	assert(loci);
	Chromosome recombinant(loci);
	
	Chromosome * currChrom = this;
	Chromosome * otherChrom = &secondChrom;
	Chromosome * tempChrom;
	
	// randomly choose a chromosome to start copying from
	if(generator->drand() < 0.5){
		currChrom = &secondChrom;
		otherChrom = this;
  	}
  
	// at each locus check for crossover event and
	// switch which chromosome is being used to copy
	// the alleles
	unsigned int totalLoci = LociCount();
	unsigned int lastLocus = totalLoci-1;
	for(unsigned int currLoc=0; currLoc<lastLocus; currLoc++){
	    recombinant[currLoc] = (*currChrom)[currLoc];
    	if(generator->drand() < (*loci)[currLoc+1].RecombinationFraction()){
			tempChrom = currChrom;
			currChrom = otherChrom;
			otherChrom = tempChrom;
    	}  
  	}
	// set final locus, no need to check for crossover
	recombinant[lastLocus] = (*currChrom)[lastLocus];
	return recombinant;
}

// Use: Output chromsome using overloaded operator
// Arg: os -- output stream
//      chrom -- chromosome to output
// Ret: output stream
std::ostream & operator << (std::ostream & os, Chromosome & chrom){
	unsigned int numLoci = chrom.LociCount();
	for(unsigned int i=0; i<numLoci; i++)
    	os << chrom[i]+1 << " ";
  	os << std::endl;
  	return os;
}



void Chromosome::WritePedFormat(std::ostream& os, uint first, uint last) {
	//std::string vals[] = { "1 1 ", "1 2 ", "2 2 " };
	unsigned int numLoci = LociCount();
	if (last == 0 || last > numLoci)
		last = numLoci;

	for(unsigned int i=first; i<last; i++)
    	os << chrom[i]+1<< " ";
  	os << std::endl;
}



void Chromosome::SetValue(uint locus, uint value) {
	assert(chrom.size() > locus);
	chrom[locus] = (value==2);
}

}
