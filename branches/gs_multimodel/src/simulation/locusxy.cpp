//
// C++ Implementation: locusxy
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locusxy.h"
#include <iomanip>
#include "utility/strings.h"

using namespace std;

namespace Simulation {

/**
 * @brief streaming operator
 */
std::ostream &operator<<(std::ostream &os, LocusXY &loc) {
	os<<loc.label<<"\t"
		<<loc.alleleFreq1<<"\t"			///<This is the allele frequency for Allele 1 X
		<<loc.alleleFreq2<<"\t"			///<This is the allele frequency for Allele 1 Y
		<<loc.mapDistance<<"\t"			///<Map distance for the X
		<<loc.mapDistanceY<<"\t"		///<Map Distance for the Y
		<<loc.location<<"\t"			///<Location for the X
		<<loc.locationY<<"\t"			///<Location for the Y
		<<loc.desc;						///<This is used to determine which locus type this is
	return os;

}
std::istream &operator>>(std::istream &is, LocusXY &loc) {
	is>>loc.label
		>>loc.alleleFreq1			///<This is the allele frequency for Allele 1 X
		>>loc.alleleFreq2			///<This is the allele frequency for Allele 1 Y
		>>loc.mapDistance			///<Map distance for the X
		>>loc.mapDistanceY		///<Map Distance for the Y
		>>loc.location			///<Location for the X
		>>loc.locationY			///<Location for the Y
		>>loc.desc;						///<This is used to determine which locus type this is
	string type = Utility::ToUpper(loc.desc.c_str());
	loc.SetType(type.c_str());
	return is;
}


void LocusXY::SetType(const char *type) {
	
	if (strcmp(type, "PAR")==0) 
		this->type = PAR;
	else if (strcmp(type, "X_ONLY")==0)
		this->type = X_Only;
	else if (strcmp(type, "Y_ONLY")==0)
		this->type = Y_Only;
	else if (strcmp(type, "X_HOMOLOG")==0)
		this->type = X_Homolog;
	else if (strcmp(type, "Y_HOMOLOG")==0)
		this->type = Y_Homolog;
	else
		return;
	desc = type;
}


void LocusXY::WriteFormatted(std::ostream &os, uint width) {
	//os<<setw(width)<<id;
	os<<setw(width-1)<<left<<label<<" ";
	os<<setw(width-1)<<left<<setprecision(6)<<alleleFreq1<<" ";
	os<<setw(width-1)<<left<<setprecision(6)<<alleleFreq2<<" ";
	os<<setw(width-1)<<left<<setprecision(width-5)<<mapDistance<<" ";
	os<<setw(width-1)<<left<<setprecision(width-5)<<mapDistanceY<<" ";
	os<<setw(width-1)<<left<<location<<" ";
	os<<setw(width-1)<<left<<locationY<<" ";
	os<<setw(width-1)<<left<<desc<<"\n";
}

void LocusXY::WriteHeader(std::ostream& os, uint width) {
	os<<setw(width-1)<<left<<"Label"<<" ";
	os<<setw(width-1)<<left<<"Freq Al1 X"<<" ";
	os<<setw(width-1)<<left<<"Freq Al1 Y"<<" ";
	os<<setw(width-1)<<left<<"Map Dist. X"<<" ";
	os<<setw(width-1)<<left<<"Map Dist. Y"<<" ";
	os<<setw(width-1)<<left<<"Position X"<<" ";
	os<<setw(width-1)<<left<<"Position Y"<<" ";
	os<<setw(width-1)<<left<<"Description"<<"\n";
}
float LocusXY::GetMapPositionX() {
	return mapPosition;
}

float LocusXY::GetMapPositionY() {
	return mapLocationY;
}
void LocusXY::SetMapPositionX(float pos) {
	mapPosition = pos;
}

int LocusXY::GetInitialValue(Utility::Random& rnd, bool isY) {
	if (isY) {
		if (type == X_Only || type == X_Homolog)
			return fixedValue;
		else
			return rnd.drand() > alleleFreq2;		
	} else {
		if (type == Y_Only || type == Y_Homolog)
			return fixedValue;
		else
			return rnd.drand() > alleleFreq1;
	}
}
	
void LocusXY::SetMapPositionY(float pos) {
	mapLocationY = pos;
}

void LocusXY::SetMapDistanceX(float dist) {
	mapDistance = dist;
}

void LocusXY::SetMapDistanceY(float dist) {
	mapDistanceY = dist;
}

float LocusXY::GetMapDistanceX() {
	return mapDistance;
}

float LocusXY::GetMapDistanceY() {
	return mapDistanceY;
}

void LocusXY::SetLocationX(int loc) {
	location = loc;
}

void LocusXY::SetLocationY(int loc) {
	locationY = loc;
}

int LocusXY::GetLocationX() {
	return location;
}

int LocusXY::GetLocationY() {
	return locationY;
}

float LocusXY::Freq1X() const {
	return alleleFreq1;
}

float LocusXY::Freq2X() const {
	return 1.0 - alleleFreq1;
}

float LocusXY::Freq1Y() const {
	return alleleFreq2;
}

float LocusXY::Freq2Y() const {
	return 1.0 - alleleFreq2;
}
bool LocusXY::operator <( const Locus& other) const {
	return location < other.GetLocation();
}
/**
 * @brief return the minor allele freq for X
 */
float LocusXY::GetMinAlleleFreqX() const {
	if (alleleFreq1 < 0.5)
		return alleleFreq1;
	else 	
		return 1.0 - alleleFreq1;
}

/**
 * @brief return the minor allele freq for Y
 */
float LocusXY::GetMinAlleleFreqY() const {
	if (alleleFreq2 < 0.5)
		return alleleFreq2;
	else 	
		return 1.0 - alleleFreq2;
}

bool LocusXY::Distill(Locus& other) {
	other = Locus(chromID, id - 1);
	other.SetLabel(label.c_str());
	if (type == Y_Homolog || type == Y_Only) {
		other.AssignFreq(alleleFreq2, 1.0 - alleleFreq2);
		other.SetLocation(GetLocationY());
		other.SetMapDistance(GetMapDistanceY());
		other.MapPosition(GetMapPositionY());
	} else if (type == X_Homolog || type == X_Only) {
		other.AssignFreq(alleleFreq1, 1.0 - alleleFreq1);
		other.SetLocation(GetLocationX());
		other.SetMapDistance(GetMapDistanceX());
		other.MapPosition(GetMapPositionX());
	} else {
		float avgAlFreq = alleleFreq1 * 0.75 + alleleFreq2 * 0.25;
		other.AssignFreq(avgAlFreq, 1.0 - avgAlFreq);
		other.SetLocation(GetLocationX());
		other.SetMapDistance(GetMapDistanceX());
		other.MapPosition(GetMapPositionX());
	}
	return true;
}

bool LocusXY::Clone(Locus& other, int t, int &fixedCount) {
	other = Locus(chromID, id);
	if (t == 0 && (type == Y_Homolog || type == Y_Only) )
		return false;
	if (t == 1 && (type == X_Homolog || type == X_Only) )
		return false;
	if (t > 1 || t < 0)
		return false;
	if ((t == 0) && (Freq1X() == 1.0 || Freq1X() == 0.0)) {
		fixedCount++;
		return false;
	}
	if (t == 1 && (Freq1Y() == 1.0 || Freq1Y() == 0.0)) {
		fixedCount++;
		return false;
	}

	other.SetLabel(label.c_str());
	//X
	if (t == 0) {
		other.AssignFreq(alleleFreq1, 1.0 - alleleFreq1);
		other.SetLocation(GetLocationX());
		other.SetMapDistance(GetMapDistanceX());
		other.MapPosition(GetMapPositionX());
	}
	else {
		other.AssignFreq(alleleFreq2, 1.0 - alleleFreq2);
		other.SetLocation(GetLocationY());
		other.SetMapDistance(GetMapDistanceY());
		other.MapPosition(GetMapPositionY());
	}
//Some properties are unaccessible due to the template-waiting to see if they are necessary, since
//this copy is really just used for LD
//	other.desc			= desc;
//	other.idx			= idx;
//	other.location		= location;
//	other.mapPosition	= mapPosition;

	return true;
}


bool LocusXY::PassThresholdMAF(float threshold) {
	return GetMinAlleleFreqX() > (threshold- 0.0001) && GetMinAlleleFreqY() > (threshold- 0.0001);
}
}
