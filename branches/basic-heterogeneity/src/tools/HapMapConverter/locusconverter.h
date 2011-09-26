//
// C++ Interface: locusconverter
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TOOLSLOCUSCONVERTER_H
#define TOOLSLOCUSCONVERTER_H

namespace Tools {

/**
Converts a given locus 0,1 value to the appropriate letter

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class LocusConverter{
public:
    LocusConverter();
    ~LocusConverter();

	string label;
	uint position;
	char *alleles;
	bool doWriteToFile;

	bool Parse(istream& file);
	string GetGenotype(int al1, int al2);
};

inline
LocusConverter::LocusConverter() : label(""),position(0), doWriteToFile(false) {
	alleles = new char[3];
}

inline
LocusConverter::~LocusConverter() { 
	delete[] alleles;
}

inline
string LocusConverter::GetGenotype(int al1, int al2) {
	assert(al1 < 2 && al2 < 2);
	stringstream ss;
	ss<<alleles[al1]<<alleles[al2];
	return ss.str();
}

inline
bool LocusConverter::Parse(istream& ss) {
	bool success=false;
	string pos;
	if (!ss.eof())	{
		string lbl;
		ss>>lbl>>pos>>alleles[0]>>alleles[1];
		label=lbl;
		position = atoi(pos.c_str());
		success=position>0;
	}
	return success;
}

}

#endif
