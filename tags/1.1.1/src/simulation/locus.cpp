// Locus.cpp

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

#include <sstream>
#include <iomanip>
#include "locus.h"

namespace Simulation {

using namespace std;

std::ostream & operator << (std::ostream & os, Locus &loc) {
	//os<<loc.id<<"\t"
	os<<loc.label<<"\t"<<loc.alleleFreq1<<"\t"<<loc.alleleFreq2<<"\t"<<loc.mapDistance<<"\t"<<loc.location<<"\t"<<loc.desc;
	return os;
}

std::istream &operator>>(std::istream &is, Locus &loc) {
	//We don't want to read past the eoln
	char line[4096];
	is.getline(line, 4096);			///<The chromosome id
	stringstream data(line);

	//data>>loc.id;
	data>>loc.label>>loc.alleleFreq1>>loc.alleleFreq2;
	data>>loc.mapDistance>>loc.location>>loc.desc;
	return is;
}

uint Locus::GetID() {
	return id;
}

uint Locus::GetChromID() {
	return chromID;
}


void Locus::SetID(uint id) {
	this->id = id;
}



void Locus::SetLabel(int chromID, int id) {
	stringstream st;
	st<<"RL"<<chromID<<"-"<<id;
	label = st.str();
}

void Locus::SetLabel(const char *label) {
	this->label = label;
}

string Locus::GetLabel() const{
	return label;
}

float Locus::GetMinAlleleFreq() const {
	float maf = Freq1();
	
//	cout<<"Evaluating MAF("<<label<<") = "<<Freq1()<<", "<<Freq2();
	
	if (Freq2() < maf)
		maf = Freq2();

//	cout<<" = "<<maf<<"\n";
	return maf;
}

bool Locus::operator<=(const Locus& other) const {
	return mapPosition <= other.mapPosition;
}

bool Locus::operator>=(const Locus& other) const {
	return mapPosition >= other.mapPosition;
}

bool Locus::operator>(const Locus& other) const {
	return mapPosition > other.mapPosition;
}

bool Locus::operator <( const Locus& other) const {
	return mapPosition < other.mapPosition;
}

void Locus::WriteHeader(std::ostream &os, uint width) {
	//os<<setw(width)<<"Id";
	os<<setw(width-1)<<left<<"Label"<<" ";
	os<<setw(width-1)<<left<<"Freq Al1"<<" ";
	os<<setw(width-1)<<left<<"Freq Al2"<<" ";
	os<<setw(width-1)<<left<<"Map Dist."<<" ";
	os<<setw(width-1)<<left<<"Position"<<" ";
	os<<setw(width-1)<<left<<"Description"<<"\n";
	
}

void Locus::WriteFormatted(std::ostream &os, uint width) {
	//os<<setw(width)<<id;
	os<<setw(width-1)<<left<<label<<" ";
	os<<setw(width-1)<<left<<setprecision(6)<<alleleFreq1<<" ";
	os<<setw(width-1)<<left<<setprecision(6)<<alleleFreq2<<" ";
	os<<setw(width-1)<<left<<setprecision(width-5)<<mapDistance<<" ";
	os<<setw(width-1)<<left<<location<<" ";
	os<<setw(width-1)<<left<<desc<<"\n";

}

void Locus::WriteMendelFormat(ostream& os) {
	os<<label<<",\t"<<chromID<<",\t"<<location<<"\n";
}
void Locus::WriteMarkerInfo(ostream &os, uint offset /*=0*/) {
	os<<label<<"\t"<<location<<"\t"<<desc<<"\n";
}

float Locus::Freq1() const {
	return alleleFreq1;
}

void Locus::SetFreq1(const float &f) {
	alleleFreq1=f;
}

float Locus::Freq2() const {
	return alleleFreq2;
}

void Locus::SetFreq2(const float &f) {
	alleleFreq2=f;
}

uint Locus::GetLocation() const {
	return location;
}

void Locus::SetLocation(uint loc) {
	location=loc;
}


bool Locus::Valid() {
	return (alleleFreq1 + alleleFreq2) > 0.9 && (alleleFreq1 + alleleFreq2) < 1.00001 && mapDistance > 0.0;
}

string Locus::GetDescription() {
	return desc;
}

void Locus::SetDescription(const char *d){
	desc = d;
}

}
