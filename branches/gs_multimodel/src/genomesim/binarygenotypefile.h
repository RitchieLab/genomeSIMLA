//
// C++ Interface: binarygenotypefile
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIMBINARYGENOTYPEFILE_H
#define GENOMESIMBINARYGENOTYPEFILE_H

#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <fstream>

namespace GenomeSIM {

namespace Formatters {

using namespace std;

/**
Parses binary genotype files

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/

#define CURR_VERSION 1
class GenotypeFile {
public:
	GenotypeFile()  { }
	virtual ~GenotypeFile() { }
	
	/**
	 * @brief Opens file for reading
	 * @param filename Name of file
	 * @param fileDetails The details regarding the file will be written here
	 */
	virtual bool Open(const char *filename) = 0;

	/**
	 * @brief closes down the filedetails	
	 */
	virtual bool Close() = 0;
};

class AsciiGenotypeFileWriter : GenotypeFile {
public:
	bool Open(const char *filename);
	bool Close();
};

class BinaryGenotypeFile : GenotypeFile {
public:

struct BinaryIndividual {
	uint individualID;						///<Unique number identifying the individual
	uint status;							///<Status
	int genotypeCount;						///<the number of SNPs 

	/**
	 * @brief Initialize so the memory is allocated locally
	 * @note If the behavior should change, you must create a new instance!
	 */
	BinaryIndividual(uint genotypeCount) ;

	~BinaryIndividual() 						 { 	}

	/**
	 * @brief Copy genotype values to the array
	 * @param genotypes Array of integers properly sized and allocated by caller
	 */
	void GetGenotypes(int *genotypes);
	int GetGenotype(uint idx);
	void SetGenotype(uint idx, int genotype);
	void SetStatus(bool isAffected=true);

	bool Load(istream &os);
	bool Save(ostream &os);

	static int AffectedValue;
	static int UnaffectedValue;

	static int Value_Missing;
	static int Value_AA;
	static int Value_Aa;
	static int Value_aA;
	static int Value_aa;
protected:
	boost::dynamic_bitset<> locusData;


	uint bitsPerBlock;
	uint totalBits;	
	uint totalBlocks;						///<

};

struct Details {
	uint version;							///<The version of the file that has been opened
	uint type;								///<The type (Case/Control or pedigree...maybe others)
	uint snpcount;							///<Number of snps (columns)
	uint indcount;							///<Number of individuals

	
	Details() : version(CURR_VERSION), type(0), snpcount(0), indcount(0) { }
	Details(uint v, uint t, uint sc, uint ic) : version(v), type(t), snpcount(sc), indcount(ic) {}
	~Details() {	 }

	bool Load(istream &os) {
		os.read((char*)&version, 4);
		os.read((char*)&type, 4);
		os.read((char*)&snpcount, 4);
		os.read((char*)&indcount, 4);	
		return true;	
	}
	bool Save(ostream &os) {
		os.write((char*)&version, 4);
		os.write((char*)&type, 4);
		os.write((char*)&snpcount, 4);
		os.write((char*)&indcount, 4);
		return true;
	}
};
	uint CurrentVersion;							///<This represents the current version of the class

    BinaryGenotypeFile(Details& details); 
    virtual ~BinaryGenotypeFile() { 	}

	/**
	 * @brief Opens file for reading
	 * @param filename Name of file
	 * @param fileDetails The details regarding the file will be written here
	 */
	virtual bool Open(const char *filename) = 0;


	/**
	 * @brief closes down the filedetails	
	 */
	bool Close() {
		return true;
	}
protected:
	void ReadDetails() {
//		uint v = 0,t = 0,sc = 0,ic =0;


		filedetails.Load(file);
		firstIndividual = currentIndividual = file.tellg();
	}

	
	void Reset() {
		if (file) 
			file.seekg(firstIndividual);
	}

	Details &filedetails;								///<Details regarding the opened file (if there is one)
	std::fstream file;										///<This is the actual file 
	std::streampos firstIndividual;							///<For resetting the file position
	std::streampos currentIndividual;						///<The current individual being read

};

class BinaryGenotypeFileWriter : BinaryGenotypeFile {
	BinaryGenotypeFileWriter(Details& details) : BinaryGenotypeFile(details) { }
	/**
	 * @brief Opens file for reading
	 * @param filename Name of file
	 * @param fileDetails The details regarding the file to be written
	 */
	bool Open(const char *filename);

	/**
	 * @brief Writes an individual to the file
	 * @return False if the file has been satisfied (based on the fileDetails)
	 */
	bool AddIndividual(BinaryGenotypeFile::BinaryIndividual* ind);
};

class BinaryGenotypeFileReader : BinaryGenotypeFile {
public:
	BinaryGenotypeFileReader(Details& details) : BinaryGenotypeFile(details) {}

	/**
	 * @brief Opens file for reading
	 * @param filename Name of file
	 * @param fileDetails The details regarding the file will be written here
	 */
	bool Open(const char *filename);

	/**
	 * @brief Resets the first individual found in the file
	 * @note Returns NULL if the file is not opened or there are no entries found
	 */
	BinaryIndividual *GetFirst();

	/**
	 * @brief Returns the first individual of the file
	 * @note Returns true if there were no errors encountered
	 */
	bool GetFirst(BinaryIndividual& i);

	/**
	 * @brief Returns the next individual
	 * @note Returns true if there were no errors encountered
	 */
	bool GetNext(BinaryIndividual& i);

};

inline
BinaryGenotypeFile::BinaryIndividual::BinaryIndividual(uint genotypeCount) : 
		individualID(0), 
		status(0), 
		genotypeCount(genotypeCount) {

	//We need to know how many blocks we want to read
	//First, the number of bits / block in our bitreader
	bitsPerBlock = boost::dynamic_bitset<>::bits_per_block;
	
	//How many bits? 3 is for the number of bits we'll, since phased information 
	//is relavent
	totalBits = genotypeCount * 3;

	//We want 1 extra bit for the status
	totalBlocks = (1 + totalBits) / bitsPerBlock;
	if (totalBits % bitsPerBlock > 0)
		totalBlocks++;
}
}
}
#endif
