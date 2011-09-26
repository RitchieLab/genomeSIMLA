// Locus.h

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
// Contains information for a single locus.
// Each locus contains allele frequencies, recombination fraction,
// direction of genotype errors, and genotype error rate
//
/////////////////////////////////////////////////////////////////////

#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <iostream>
#include <assert.h>
#include "utility/random.h"
#include <math.h>

namespace Simulation {

class Locus{
  public:
    // Constructors
	/**
	 * @brief Constructor
	 * @param errorDirect The direction any error will occur
	 * @param recombFreq The likelihood that we'll jump chromosomes at this locus during crossing
	 * @param alFreq1 Frequency of allele 1
	 * @param alFreq2 Frequency of allele 2
	 * @param errRate The rate at which an error will occur at this locus
	 */
    Locus(float recombFrac, float errRate, uint id, uint location) : 
			recombFraction(recombFrac), errorRate(errRate), location(location), desc("") {
		SetLabel(id);
	}

	/**
	 * This is to be used only when loading from file
	 */
	Locus() : alleleFreq1(0.0), alleleFreq2(0.0), recombFraction(0.0), errorRate(0.0), location(0), desc("") {
		SetLabel(0);
	}


	void SetLabel(int id);
	void SetLabel(const char *label);
	std::string GetLabel();

	/**
	 * @brief Get/set error rate
	 */
    float ErrorRate() const{return errorRate;}
	void ErrorRate(float err) { errorRate=err; }
    
	/**	
 	 * @brief return the frequency at allele, allele. 
	 * @param allele indicate which allele you are referring to (1 or 2)
	 */
	float GetAlleleFreq(uint allele);    
	
	/**
	 * @brief Returns the local recombination fraction
	 */
    double RecombinationFraction(){return recombFraction;}

	/**
	 * @brief Setup the allele frequencies randomly
	 */
	void RandomizeFreq(float min, float max);

	/**
	 * @brief Assign Frequencies exactly according to the parameters passed
	 */
	void AssignFreq(float al1, float al2);

	/**
	 * @brief Returns the frequency for the first allele
	 */
	float Freq1();

	/**
	 * @brief Returns the frequency for the second allele 
 	 */
	float Freq2();
	
	/**
	 * @brief streaming operator
	 */
    friend std::ostream &operator<<(std::ostream &os, Locus &loc);
	friend std::istream &operator>>(std::istream &is, Locus &loc);

	/**
	 * @brief Write the header details to the stream
	 * @param os the stream
	 * @Param width the amount of space between the columns
	 */
	void WriteHeader(std::ostream &os, uint width);

	/**
	 * @brief Write the information of the locus to the stream
	 * @param os the stream
	 * @Param width the amount of space between the columns
	 */
	void WriteFormatted(std::ostream &os, uint width);

	/**
	 * @brief return the physical location of the locus on the chromosome
	 */
	uint GetLocation();
		
	/**
	 * @brief set the physical location of the locus on the chromosome
	 */
	void SetLocation(uint loc);

	/**
	 * @brief Verifies that the locus has reasonable frequency / fraction information 
	 */
	bool Valid();
	
	
	/**
	 * @brief Sets the index of the snp on the chromosome
	 */
	//void SetID(uint id);

	/**
	 * @brief Write the marker information to the stream	
 	 */
	void WriteMarkerInfo(std::ostream &os);

	/**
	 * @brief Return the number of base pairs between this node and the one before it
	 */
	uint GetBasePairDistance();

	/**
	 * @brief return the description	 
	 */
	std::string GetDescription();
	void SetDescription(const char *description);
private:
	float alleleFreq1;						///<First allele
	float alleleFreq2;						///<Second allele

    float recombFraction,					///<Recombination fraction
		errorRate;							///<Rate at which errors occur

	//int id;									///<The snp ID for this location
	std::string label;
	int location;							///<The base-pair location of the locus on 
	std::string desc;						///<Allow for something to annotate this locus
};



inline
float Locus::GetAlleleFreq(uint allele) {
	float af = 0.0;
	if (allele == 1)
		af = alleleFreq1;
	else if (allele == 2)	
		af = alleleFreq2;
	else {
		std::cout<<"Allele index should be one based\n";
		abort();
	}
	return af;
}



inline
uint Locus::GetBasePairDistance() {
	//By default, we will calculate kosambi.
	double Cab2 = 2*recombFraction;
	double cmCount=0.25 * log((1.0 + Cab2)/(1.0-Cab2));
	return (uint)(cmCount * 1000000.0);
}

inline
void Locus::AssignFreq(float al1, float al2) {
	alleleFreq1 = al1;
	alleleFreq2 = al2;
}

inline
void Locus::RandomizeFreq(float min, float max) {
	assert(min <= max);
	float diff = max-min;
	float value = min + Utility::Random::globalGenerator.drand() * diff;
	
	if (min > 0.5) {
		alleleFreq2 = 1.0 - value;
		alleleFreq1 = value;
	}
	else {
		alleleFreq1 = 1.0 - value;
		alleleFreq2 = value;
	}

}

}




#endif 
