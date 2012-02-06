//
// C++ Implementation: allelicdecoder
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "allelicdecoder.h"
#include "utility/strings.h"
#include "utility/random.h"
#include <sstream>
#include <iostream>

namespace Genetics {
void AllelicDecoder::LocusParser::Seed(string &data) {
	locus[data]++;
	seedCount++;
}
int AllelicDecoder::LocusParser::Decode(string &encoding) {
	int value = missingValue;
	assert(atoi(encoding.c_str()) < 10);

	if (encoding == missingData)
		return InferValue();

	//Here, we'll assign 0 to missing and 1 to allele one, 2 to allele two etc...
	for (uint i=0; i<maxPoly && value == missingValue; i++) 
		if (encodings[i] == encoding)
			value = i + 1;

		
	//Eventually, we can return an index into the real decoded values
	return value;
}

int AllelicDecoder::LocusParser::InferValue() {
	//For now, we are assuming biallelic 
	//WE have to update this if we really do encounter more alleles
	assert(locus.size() <= 3);
	if (Random::globalGenerator((int)seedCount) < locus[encodings[0]])
		return 1;
	else
		return 2;
}

bool AllelicDecoder::LocusParser::Verify(uint alleleCount) {
	map<int, string> invertedMap;
	CodeMap::iterator itr = locus.begin();
	CodeMap::iterator end = locus.end();

	while (itr != end) {
		invertedMap[itr->second] = itr->first;
		itr++;
	}

	map<int, string>::reverse_iterator i = invertedMap.rbegin();
	map<int, string>::reverse_iterator e = invertedMap.rend();
	//Because invertedMap is sorted by key (which happens to be counts) we find that we are marching down the list
	encodings.clear();
	while (i!=e) {
		string key = i->second;
		//We don't want to mess with missing data
//		if (key != missingData) {
			encodings.push_back(key);
//		}
		i++;

	}
	maxPoly = alleleCount;	
//	assert(locus.size() <= alleleCount);
	return locus[missingData] == 0;
}

/**
 * @brief Feed all possible lines to the decoder so it can determine which is the major allele.
 * @param line This contains all encoded alleles for at least 1 chromosome
 * @param stride The number of columns per genotype (diploid on 1 line would be 2)
 */
void AllelicDecoder::SeedEncodings(const char *line, int stride, int headerColumns) {
	int columns = CountColumns(line);
	if (columns < 1)
		return;

	assert(columns > headerColumns);
	columns = (columns - headerColumns) / stride;

	if (loci.size() < (uint)columns)
		loci.resize(columns, LocusParser(missingData.c_str()));	

	string temp;
	stringstream ss(line);
	//Skip over headers
	int n=0;
	while (n++ < headerColumns)
		ss>>temp;

	for (int i=0; i<columns; i++) {
		for (int j=0; j<stride; j++) {
			ss>>temp;
			loci[i].Seed(temp);
		}
	}
	
}

/**
 * @brief Ensure that we don't have any errors, and we recognize missing data properly
 * @param maxPoly The maximum number of variant alleles
 * @param missingEncoding The item that denotes missing data
 */
bool AllelicDecoder::Verify(uint maxPoly) {
	uint count = loci.size();
	bool success = true;
	hasMissingLoci.clear();
	hasMissingLoci.resize(count, false);
	for (uint i=0; i<count && success; i++) {
		hasMissingLoci[i] = loci[i].Verify(maxPoly);
	}
	cout<<"Verification of seed data has completed. Total SNPs: "<<count<<" #with Missing Data: "<<hasMissingLoci.count()<<"\n";
	return success;
}


int AllelicDecoder::GetMissingValue() {
	return -1;
}

/**
 * @brief Actually translate an array of characters to a single chromosome
 * @param line character encoded allelic data
 * @param stride number of columns per snp
 * @param which column (above) we are reading from 
 * @note This just lets us skip over lines where both paternal and maternal are on same line
 */
bool AllelicDecoder::ParseLine(const char *line, vector<int>& chromosome, int stride, int headerCols, int offset) {
	stringstream ss(line);
	string temp;

	int pos = 0;
	int lociCount = loci.size();

	//We need to make sure nothing is left behind from the last time
	chromosome.clear();
	chromosome.resize(lociCount, 0);

	while (headerCols-- > 0)
		ss>>temp;

	for (int i=offset; i<lociCount;) {
		ss>>temp;
		
//		if (!hasMissingLoci[i]) 
			chromosome[pos] = loci[pos].Decode(temp);
//		else
//			chromosome[pos] = 1;
		pos++;
	
		i+=stride;
	}
	return lociCount > offset;
}
AllelicDecoder::AllelicDecoder::AllelicDecoder(const char *missingData) : missingData(missingData)
{
}


AllelicDecoder::~AllelicDecoder()
{
}


}
