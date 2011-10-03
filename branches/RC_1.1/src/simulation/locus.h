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
    Locus(float recombFrac, uint chromID, uint id, uint location) : 
			id(id), chromID(chromID) ,alleleFreq1(0.0), alleleFreq2(0.0), mapDistance(recombFrac), location(location), desc(""), mapPosition(0.0) {
		SetLabel(chromID, id);
	}

	/**
	 * This is to be used only when loading from file
	 */
	Locus(uint chromID, uint id) : id(id), chromID(chromID), alleleFreq1(0.0), alleleFreq2(0.0), mapDistance(0.0), location(0), desc(""), mapPosition(0.0) {
		SetLabel(chromID, id);
	}
		
	//This is to be used by STL only or when we are about to assign it from something else
	Locus() :id(0), chromID(0) ,alleleFreq1(0.0), alleleFreq2(0.0), mapDistance(0.0), location(0), desc(""), mapPosition(0.0) {}

	virtual ~Locus() {
		//std::cout<<"*\tLocus: "<<id<<" is going away\n";
	}


	void SetLabel(int chromID, int id);
	void SetLabel(const char *label);
	std::string GetLabel() const;

	/**
	 * @brief Get/set error rate
	 */
//    float ErrorRate() const{return errorRate;}
//	void ErrorRate(float err) { errorRate=err; }
    
	/**	
 	 * @brief return the frequency at allele, allele. 
	 * @param allele indicate which allele you are referring to (1 or 2)
	 */
	float GetAlleleFreq(uint allele);
	
	/**
	 * @brief Returns the local recombination fraction
	 */
    double MapDistance(){return mapDistance;}
	double MapPosition() {return mapPosition; }
	void MapPosition(double pos) { mapPosition=pos; }

	void SetMapDistance(const double &dist) { mapDistance = dist; }

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
	float Freq1() const;

	void SetFreq1(const float &f);

	/**
	 * @brief Returns the frequency for the second allele 
 	 */
	float Freq2() const;
	
	void SetFreq2(const float &f);

	int GetInitialValue(Utility::Random& rnd, bool isY);

	/**
	 * @brief create a pseudo copy of the local locus into loc. This is useful for downgrading a derived class
	 * @param loc Object to be written to
	 * @param t type of copy. For XY, this can be used to write either X or Y loci
	 * @param Return True if the locus of type, t, makes sense
	 */
	virtual bool Clone(Locus *loc, int t, int &fixedCount);

	/**
	 * @Brief Similar to clone, except in some cases, it might blend certain features (XY allele frequencies)
	 */
	virtual bool Distill(Locus& loc);
	/**
	 * @brief streaming operator
	 */
    friend std::ostream &operator<<(std::ostream &os, Locus &loc);
	friend std::istream &operator>>(std::istream &is, Locus &loc);

	bool operator<(const Locus& other) const;
	bool operator>(const Locus& other) const;
	bool operator<=(const Locus& other) const;
	bool operator>=(const Locus& other) const;

	/**
	 * @brief Write the header details to the stream
	 * @param os the stream
	 * @Param width the amount of space between the columns
	 */
	virtual void WriteHeader(std::ostream &os, uint width);

	/**
	 * @brief Write the information of the locus to the stream
	 * @param os the stream
	 * @Param width the amount of space between the columns
	 */
	virtual void WriteFormatted(std::ostream &os, uint width);

	/**
	 * @brief return the physical location of the locus on the chromosome
	 */
	uint GetLocation() const;
		
	/**
	 * @brief set the physical location of the locus on the chromosome
	 */
	void SetLocation(uint loc);

	/**
	 * @brief Verifies that the locus has reasonable frequency / fraction information 
	 */
	virtual bool Valid();
		
	/**
	 * @brief Sets the index of the snp on the chromosome
	 */
	uint GetID();

	uint GetChromID();

	void SetID(uint id);
	/**
	 * @brief Write the marker information to the stream
	 * @param  os Stream to be written to
	 * @param offset the position for the first marker of the chromosome 
  	 * 				 (if not already figured correctly for the chromosome)
 	 */
	void WriteMarkerInfo(std::ostream &os, uint offset = 0);
	void WriteMendelFormat(std::ostream& os);
	/**
	 * @brief Return the number of base pairs between this node and the one before it
	 */
	uint GetBasePairDistance();

	/**
	 * @brief return the description	 
	 */
	std::string GetDescription();
	void SetDescription(const char *description);


	float GetMinAlleleFreq() const;

	/**
	 * @Brief allow the object to determine if it's MAF passes the threshold (is > thresh)
	 */
	bool PassThresholdMAF(float threshold);

	virtual bool operator()(const Locus  s2) {
		if (location == s2.location)
			return idx < s2.idx;
		else
			return location < s2.location;
	}	

	Locus(const Locus& other) : idx(other.idx), id(other.id), chromID(other.chromID), 
				alleleFreq1(other.alleleFreq1), alleleFreq2(other.alleleFreq2), 
				mapDistance(other.mapDistance), 
				label(other.label), location(other.location), desc(other.desc), mapPosition(other.mapPosition) { }


struct PositionCalibration {
	float map;	
	int idx;
	void Calibrate(Locus& loc) {
		map+=loc.MapDistance();		
		loc.MapPosition(map);
		loc.SetID(idx++);
	}
	PositionCalibration() : map(0.0), idx(0) { }
};
protected:
	uint idx;								///<Keep track with the placement within the chromosome (not positional, but in relation to the bits)
	uint id;								///<The snp ID for this location
	uint chromID;							///<The chromosome this belongs on
	float alleleFreq1;						///<First allele
	float alleleFreq2;						///<Second allele

    float mapDistance;						///<Map Distance
//		errorRate;							///<Rate at which errors occur

	std::string label;
	int location;							///<The base-pair location of the locus on 
	std::string desc;						///<Allow for something to annotate this locus
	double mapPosition;						///<The place on the chromosome in terms of genetic distance
};

inline
bool Locus::Distill(Locus& other) {
	other = *this;
	return true;
}
inline
bool Locus::Clone(Locus* other, int t, int &fixedCount) {
	if (other)		
		delete other;
	if (alleleFreq1 < 0.00001 || alleleFreq2 < 0.00001) {
		fixedCount++;
		return false;
	}
	other = new Locus(*this);
	
	return true;
}

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
int Locus::GetInitialValue(Utility::Random& rnd, bool isY) {
	return rnd.drand() > alleleFreq1;
}

inline
bool Locus::PassThresholdMAF(float threshold) {
	return GetMinAlleleFreq() > (threshold - 0.0001);
}

inline
uint Locus::GetBasePairDistance() {
	//Map distance is in centiMorgans
	return (uint)(mapDistance * 1000000.0);
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
	float value = min + (float)Utility::Random::globalGenerator.drand() * diff;
	
	if (min > 0.5) {
		alleleFreq2 = (float)1.0 - value;
		alleleFreq1 = value;
	}
	else {
		alleleFreq1 = (float)1.0 - value;
		alleleFreq2 = value;
	}

}

}




#endif 
