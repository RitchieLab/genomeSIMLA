//
// C++ Implementation: binarygenotypefile
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "binarygenotypefile.h"

namespace GenomeSIM {
namespace Formatters {

using namespace std;
using namespace boost;

int BinaryGenotypeFile::BinaryIndividual::AffectedValue = 1;
int BinaryGenotypeFile::BinaryIndividual::UnaffectedValue = 0;
int BinaryGenotypeFile::BinaryIndividual::Value_AA = 0;
int BinaryGenotypeFile::BinaryIndividual::Value_Aa = 1;
int BinaryGenotypeFile::BinaryIndividual::Value_aA = 2;
int BinaryGenotypeFile::BinaryIndividual::Value_aa = 3;
int BinaryGenotypeFile::BinaryIndividual::Value_Missing = 4;


bool BinaryGenotypeFile::BinaryIndividual::Save(ostream &file) {
	file.write((char*)&individualID, 4);

	dynamic_bitset<>::block_type raw[totalBlocks];

	//Convert the bitset to raw data
	to_block_range(locusData, raw);

	file.write((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));	

	return true;
}

void BinaryGenotypeFile::BinaryIndividual::GetGenotypes(int *genotypes) {
	for (int i=0; i<genotypeCount; i++) 
		genotypes[i]=GetGenotype(i);
}

void BinaryGenotypeFile::BinaryIndividual::SetStatus(bool isAffected/*=true*/) {
	locusData[0] = isAffected;
}
void BinaryGenotypeFile::BinaryIndividual::SetGenotype(uint idx, int genotype) {
	if (genotype == Value_Missing) {
		locusData[idx++]=0;
		locusData[idx++]=0;
		locusData[idx++]=0;
		return;
	} else if (genotype == Value_AA) {
		locusData[idx++]=0;
		locusData[idx++]=0;
		locusData[idx++]=1;
		return;
	} else if (genotype == Value_Aa ) {
		locusData[idx++]=1;
		locusData[idx++]=0;
		locusData[idx++]=1;
		return;
	} else if (genotype == Value_aA) {
		locusData[idx++]=0;
		locusData[idx++]=1;
		locusData[idx++]=1;
		return;
	} else if (genotype == Value_aa) {
		locusData[idx++]=1;
		locusData[idx++]=1;
		locusData[idx++]=1;
		return;
	}
}

int BinaryGenotypeFile::BinaryIndividual::GetGenotype(uint idx) {
	int genotype = 0;
	//                     00            00        01        10        11
	int returnValues[] = {Value_Missing, Value_AA, Value_Aa, Value_aA, Value_aa};

	//The first bit is for status
	idx++;
	if (locusData[idx+2])
		genotype = 2 * locusData[idx] + locusData[idx++];

	return returnValues[genotype];		
}

bool BinaryGenotypeFile::BinaryIndividual::Load(istream &file) {
	file.read((char*)&individualID, 4);
	
	dynamic_bitset<>::block_type raw[totalBlocks];
	
	file.read((char*)&raw, (sizeof(dynamic_bitset<>::block_type)*totalBlocks));

	locusData=dynamic_bitset<>(raw[0], raw[totalBlocks-1]);
	//Initialize the bitset from the file
	
	if (locusData[0])
		status = AffectedValue;
	else
		status = UnaffectedValue;	
	
	return true;
}

bool AsciiGenotypeFileWriter::Open(const char *filename) {

	assert(0);
//	file.open(filename, ios_base::out);
	return true;
}

bool AsciiGenotypeFileWriter::Close() {
	assert(0);
	return true;
}

/**
 * @brief Opens file for reading
 * @param filename Name of file
 * @param fileDetails The details regarding the file to be written
 */
bool BinaryGenotypeFileWriter::Open(const char *filename) {
	file.open(filename, ios_base::out);
	ReadDetails();	
	return true;
}


bool BinaryGenotypeFileWriter::AddIndividual(Individual *ind) {
	size_t snpCount = ind->GetSnpCount();
	
	BinaryIndividual basicInd(snpCount);
	
	for (size_t i=0; i<snpCount; i++) {
		basicInd.SetGenotype( ind->GetGenotype( i ));
	}

	return basicInd->Save(file);
}
/**
 * @brief Writes an individual to the file
 * @return False if the file has been satisfied (based on the fileDetails)
 */
bool BinaryGenotypeFileWriter::AddIndividual(BinaryGenotypeFile::BinaryIndividual* ind) {
	ind->Save(file);
	return true;
}

/**
 * @brief Opens file for reading
 * @param filename Name of file
 * @param fileDetails The details regarding the file will be written here
 */
bool BinaryGenotypeFileReader::Open(const char *filename) {
	file.open(filename, ios_base::in);
	ReadDetails();	
	return true;
}

/**
 * @brief Resets the first individual found in the file
 * @note Returns NULL if the file is not opened or there are no entries found
 */
BinaryGenotypeFile::BinaryIndividual *BinaryGenotypeFileReader::GetFirst() {
	BinaryGenotypeFile::BinaryIndividual *i = new BinaryGenotypeFile::BinaryIndividual(filedetails.snpcount);
	i->Load(file);
	return i;
}

/**
 * @brief Returns the first individual of the file
 * @note Returns true if there were no errors encountered
 */
bool BinaryGenotypeFileReader::GetFirst(BinaryIndividual& i) {
	Reset();
	return i.Load(file);
}

/**
 * @brief Returns the next individual
 * @note Returns true if there were no errors encountered
 */
bool BinaryGenotypeFileReader::GetNext(BinaryIndividual& i) {
	return i.Load(file);
}
}
}
