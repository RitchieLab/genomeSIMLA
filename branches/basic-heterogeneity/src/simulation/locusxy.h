//
// C++ Interface: locusxy
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONLOCUSXY_H
#define SIMULATIONLOCUSXY_H

#include "locus.h"
#include <iostream>
namespace Simulation {

/**
@brief Specialized Locus type for XY.
This class differs from Locus in that each locus is one of several types: XOnly, YOnly, XHomolog, YHomolog and PAR. 

	@Note An important point is that all default positional details are left associated with the X. This allows the standard code to apply to XO with 2 Xs without any change in any code. However, we have specialized code (PAR) for the XY cross-overs, which is based solely on the Y. The remaining Get???X functions are left for clarities sake should the locus object be used by other classes which use it as a specialized, XY locus object. 
	
	Specialized Roles:
		X_Only		These SNPS do NOT exist in the Y which means that XX can be heterozygous, but XY can only be homozygous
		Y_Only		These SNPs do NOT exist in the X, which means that XX won't have a genotype, and XY will only be homozygous
		X_HOMOLOG	These SNPs are polymorphic in X only. This only affects initialization
		Y_HOMOLOG	These SNPs are polymorphic in Y only. This affects initialization only
		PAR			These SNPs occur in both X and Y, but for Y, it is the only region of XO
	@author 
*/
class LocusXY : public Locus
{
public:
	enum LocusType {
		Normal = 0, 
		X_Only,
		Y_Only,
		X_Homolog,
		Y_Homolog,
		PAR
	};

	 /** 
	 * @brief Constructor
	 * @param errorDirect The direction any error will occur
	 * @param recombFreq The likelihood that we'll jump chromosomes at this locus during crossing
	 * @param alFreq1 Frequency of allele 1
	 * @param alFreq2 Frequency of allele 2
	 */
    LocusXY(float recombFrac, uint chromID, uint id, uint location) : 
			Locus(recombFrac, chromID, id, location), 
			type(Normal), mapDistanceY(0.0), locationY(0), mapLocationY(0.0), fixedValue(false) { }

	/**
	 * This is to be used only when loading from file
	 */
	LocusXY(uint chromID, uint id) : 
			Locus(chromID, id), 
			type(Normal), mapDistanceY(0.0), locationY(0), mapLocationY(0.0), fixedValue(false) { }
		
	//This is to be used by STL only or when we are about to assign it from something else
	LocusXY() : type(Normal), mapDistanceY(0.0), locationY(0), mapLocationY(0.0), fixedValue(false) {}

    ~LocusXY() { }


	/**
	 * @brief streaming operator
	 */
    friend std::ostream &operator<<(std::ostream &os, LocusXY &loc);
	friend std::istream &operator>>(std::istream &is, LocusXY &loc);


	/**
	 * @brief Returns an initial value based on local frequencies and a single random draw
	 */
	int GetInitialValue(Utility::Random& rnd, bool isY);
	/**
	 * @brief Returns the map distance for X chromosomes (this and the previous snp)
	 */
	float GetMapDistanceX();
	/**
	 * @brief Returns the map distance for Y chromosomes (this and the previous snp)
	 */
	float GetMapDistanceY();

	float GetMapPositionX();
	float GetMapPositionY();
	void SetMapPositionX(float pos);
	void SetMapPositionY(float pos);
	/**
	 * @brief Sets the map distance for X chromosomes (this and the previous snp)
	 */
	void SetMapDistanceX(float dist);
	/**
	 * @brief Sets the map distance for Y chromosomes (this and the previous snp)
	 */
	void SetMapDistanceY(float dist);

	/**
	 * @brief Returns the physical location (in basepairs) for the locus (X)
	 */
	int GetLocationX();
	/**
	 * @brief Returns the physical location (in basepairs) for the locus (Y)
	 */
	int GetLocationY();

	/**
	 * @brief Sets the physical location (in basepairs) for the locus (X)
	 */
	void SetLocationX(int location);
	/**
	 * @brief Sets the physical location (in basepairs) for the locus (Y)
	 */
	void SetLocationY(int location);

	
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
	 * @brief return the minor allele freq for X
	 */
	float GetMinAlleleFreqX() const;

	/**
	 * @brief return the minor allele freq for Y
	 */
	float GetMinAlleleFreqY() const;

	float Freq1X() const;
	float Freq2X() const;
	float Freq1Y() const;
	float Freq2Y() const;

	bool Valid() { return !(alleleFreq1 == 0 && alleleFreq2 ==0 && mapDistance==0 && mapDistance==0); }
	LocusType type;

	void SetType(const char *type);
	bool operator <( const Locus& other) const;
	bool PassThresholdMAF(float threshold);

	bool Clone (Locus& other, int t, int &fixedCount);
	bool Distill(Locus &other);
struct PositionCalibration {
	float mapX;
	float mapY;
	void Calibrate(LocusXY& loc) {
		mapX+=loc.GetMapDistanceX();
		mapY+=loc.GetMapDistanceY();
//std::cerr<<"Calibrate("<<loc.label<<") "<<mapX<<" "<<mapY<<"\n";
		loc.SetMapPositionX(mapX);
		loc.SetMapPositionY(mapY);
	}
	PositionCalibration() : mapX(0.0), mapY(0.0) { }
};
protected:
    float mapDistanceY;							///<Distance for Y from previous SNP
	int locationY;								///<Location (basepairs) for Y
	float mapLocationY;							///<Location in terms of genetic distance
	bool fixedValue;							///<Fixed allele value at this locus
};



}

#endif
