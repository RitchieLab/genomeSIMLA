//
// C++ Interface: allelicdecoder
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSALLELICDECODER_H
#define GENETICSALLELICDECODER_H

#include "utility/types.h"
#include <map>
#include <vector>
#include <string>

namespace Genetics {

using namespace std;
using namespace Utility;
/**
Help to parse various allelic encodings

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class AllelicDecoder {
public:
	typedef map<string, int> CodeMap;	
	struct LocusParser {
		CodeMap locus;		///<All possible encodings

		uint maxPoly;		///<Number of polymorphic alleles		
		uint seedCount;		///<Just to keep up with how many points we've observed


		int missingValue;	///<This is what we want missing to resolve to
		string missingData;

		/** When we are done, we'll order them according to frequency (descending) */
		vector<string> encodings;	

		LocusParser(const char *missingData) : maxPoly(2), seedCount(0), missingValue(-1), missingData(missingData) { }

		void Seed(string &encoding);
	
		int Decode(string &encoding);
		int InferValue();
		/**
	 	 * @brief Checks that we don't have too many different alleles
	 	 */		 
		bool Verify(uint alleleCount);

	};



	AllelicDecoder(const char *missingData);
	~AllelicDecoder();

	/**
	 * @brief Feed all possible lines to the decoder so it can determine which is the major allele.
	 * @param line This contains all encoded alleles for at least 1 chromosome
	 * @param stride The number of columns per genotype (diploid on 1 line would be 2)
	 */
	void SeedEncodings(const char *line, int stride, int headerColumns);
								
	/**
	 * @brief Ensure that we don't have any errors, and we recognize missing data properly
	 * @param maxPoly The maximum number of variant alleles
	 * @param missingEncoding The item that denotes missing data
	 */
	bool Verify(uint maxPoly);

	/**
	 * @brief Actually translate an array of characters to a single chromosome
	 * @param line character encoded allelic data
	 * @param stride number of columns per snp
	 * @param which column (above) we are reading from 
	 * @note This just lets us skip over lines where both paternal and maternal are on same line
	 */
	bool ParseLine(const char *line, vector<int>& chrom, int stride, int headerCols, int offset);

	int GetMissingValue();
protected:
	vector<LocusParser> loci;
	BitSetType hasMissingLoci;
	string missingData;
};

}

#endif
