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
	os<<loc.label<<"\t"<<"\t"<<loc.alleleFreq1<<"\t"<<loc.alleleFreq2<<"\t"<<loc.recombFraction<<"\t"<<loc.location<<"\t"<<loc.desc;
	return os;
}

std::istream &operator>>(std::istream &is, Locus &loc) {
	//We don't want to read past the eoln
	char line[4096];
	is.getline(line, 4096);			///<The chromosome id
	stringstream data(line);

	data>>loc.label>>loc.alleleFreq1>>loc.alleleFreq2;
	data>>loc.recombFraction>>loc.location>>loc.desc;
	return is;
}


void Locus::SetLabel(int id) {
	stringstream st;
	st<<"rl-"<<id;
	label = st.str();
}

void Locus::SetLabel(const char *label) {
	this->label = label;
}

string Locus::GetLabel() {
	return label;
}


void Locus::WriteHeader(std::ostream &os, uint width) {
	os<<setw(width)<<"Label";
	os<<setw(width)<<"Rate";
	os<<setw(width)<<"Freq Al1";
	os<<setw(width)<<"Freq Al2";
	os<<setw(width)<<"Recomb. Fraction";
	os<<setw(width)<<"Description"<<"\n";
	
}

void Locus::WriteFormatted(std::ostream &os, uint width) {
	os<<setw(width)<<label;
	os<<setw(width)<<setprecision(3)<<errorRate;
	os<<setw(width)<<setprecision(6)<<alleleFreq1;
	os<<setw(width)<<setprecision(6)<<alleleFreq2;
	os<<setw(width)<<setprecision(10)<<recombFraction;
	os<<setw(width)<<location;
	os<<setw(width)<<desc<<"\n";

}


void Locus::WriteMarkerInfo(ostream &os) {
	os<<label<<"\t"<<location<<"\t"<<desc<<"\n";
}

float Locus::Freq1() {
	return alleleFreq1;
}

float Locus::Freq2() {
	return alleleFreq2;
}

uint Locus::GetLocation() {
	return location;
}

void Locus::SetLocation(uint loc) {
	location=loc;
}


bool Locus::Valid() {
	return (alleleFreq1 + alleleFreq2) > 0.9 && recombFraction > 0.0;
}

string Locus::GetDescription() {
	return desc;
}

void Locus::SetDescription(const char *d){
	desc = d;
}

}
