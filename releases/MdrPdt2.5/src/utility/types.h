//
// C++ Interface: types
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TYPES_H
#define TYPES_H

//Required for DRAND48
//#include <acconfig.h>

//#include <stdint.h>
#include <vector>
#include <string>

#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <vector>


#define MAX_LINE_LENGTH 600000
#ifdef TEST_APP
#include <cppunit/extensions/HelperMacros.h>
#endif


namespace Utility {

using namespace boost;

typedef std::vector<std::string> StringArray;

typedef boost::tokenizer<char_separator<char> > strtokenizer;

typedef boost::dynamic_bitset<> BitSetType;
typedef std::vector<BitSetType> BitSetArray;

/**
 * @brief This is the storage medium for genetypes. 
 * It provides a label and bitarray
 */
struct GtStorage {
	std::string label;						///<This is used to identify a given genotype
	BitSetType individuals;				///<Which individuals in a given snp have this genotype
	GtStorage(const char *label, BitSetType individuals) : label(label), individuals(individuals) {}		
	GtStorage() : label("") {}
};

typedef std::vector<GtStorage> GenotypeArray;

//typedef uint uint32_t;

}

#endif

