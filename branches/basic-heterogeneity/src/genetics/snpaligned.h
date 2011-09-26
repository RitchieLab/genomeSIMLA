//
// C++ Interface: snpaligned
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSSNPALIGNED_H
#define GENETICSSNPALIGNED_H

#include <iostream>
#include <vector>
#include "utility/utility.h"

namespace Genetics {

using namespace boost;
using namespace std;
using namespace Utility;



/**
@brief In memory representation of a single snp. 
Each snp will have a set of genotype oriented bitsets which indicate which individuals have a given genotype for the particular snp. 
@todo There are some rules that are required for validation of the snp.
<UL><LI>In one case, if all individuals have the same phenotype, then it should be excluded
</UL>
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class SnpAligned{
friend class SnpPool;
public:


struct Label {
	uint lineNumber;
	uint label;		

	Label(uint lineNumber, uint label) : lineNumber(lineNumber), label(label) {}
	Label() {
		lineNumber = 0;
		label = 0;
	}
};
	/**
	 * @brief Crutch for checking on the contents of binary files. It's an array of ascii values associated with the 
	 * values read in.
	 */
//	char *asciiGenotypes;

	/**
	 * @brief Import genotype data from a string of encoded genotypes
	 * @param genotype The encoded string: i.e. 02211010200 where 0 is homozygote 1, 1 is heterozygote and 2 is homozygote 2
	 * @return the number of genotypes found
	 */
	int ImportSnp(uint linenumber, uint label, const char *genotype);

	/**
	 * @brief Imports genotype data from an array of integers. 
	 * @param linenumber Indicates the line number 
	 * @param label The label # associated with the entry
	 * @param values Compressed integer values indicating the various individual's genotypes for a given snp
	 * @param count The number of integers expected to be found in values[]
	 * @param framesize The number of bits used to represent a single individual
	 */
	int ImportCompressedBinSnp(uint linenumber, uint label, uint values[], const int count = 1, const int framesize = 2);

	/**
	 * @brief Import genotype data from a string of encoded genotypes of the form: BP BP
  	 * @param linenumber The line number the entry was found on
	 * @param label The numerical label associated with the entry 
	 * @param genotype The encoded string: i.e. AA AT AT TT where 0 is homozygote, both heterozygotes and then the other homozygote
	 * @return The number of genotypes found (each reprenting and individual)
	 * @note Currently, the file has no header information: Each line is assumed to be a valid list of individual genotypes for
	 * a given SNP. 
	 * <P>It is important to also note that this function expects the data to be genetic information
	 */
	int ImportBpSnp(uint linenumber, uint label, const char *genotype);

	/**
	 * @brief Returns a reference to the bitset for a given genotypeID
	 */
	BitSetType &GetGenotype(uint genotypeID);	

	/**
	 * @brief Returns the allelic data associated with a given genotype. 
	 * @note Currently, this is hardcoded into genetic snps. We need to discuss the approach and requirements for environmental encoding
 	 */
	string GetGenotypeAlleles(uint genotypeID);


	/**
	 * @brief Returns the inclusion status for the snp (such as INCLUDED or EXCLUDED-CONSTANT)
	 */
	const char *GetInclusionStatus();

	uint GetLabel(uint labelidx);

	uint GetLineNumber(uint labelidx);

	uint GetLabelCount();

	/**
	 * @brief Returns the number of individuals associated with a given genotypeID
	 */
	uint CountIndividuals(uint genotypeID);

	/**
	 * @brief Sets the individual count- required for parsing compressed binary array data
	 */
	void SetIndividualCount(int individualCount);

	/**
	 * @brief Returns the number of genotypes associated with this snp
	 */
	uint CountGenotypes();

	/**
	 * @brief Formulates the a new snp representing the punett product of the snp and the parameter
 	 */
	SnpAligned *punett(SnpAligned* other);

	/**
	 * @brief Returns the Maximum Class Difference Metric (MD) for the current SNP. 
	 * @param mask This represents the bit mask where the ones indicate the affected individuals
	 */
	//float MaxDifferential(const CaseControlStatus &mask);

	/**
	 * @brief Returns the Maximum Class Difference Metric for the select genotype
	 * @param genotype The genotype to be considered
	 * @param mask The mask used to determine which individuals are affected 
	 */
	//float MaxDifferential(uint genotype, const CaseControlStatus &mask);

	/**
	 * @brief Calculates the difference if ALL genotypes completely differentiate the affected from unaffected
	 */
	int IsPerfect(const BitSetType &mask);

	/**
	 * @brief Returns the perfect Difference only if all individuals are differentiated according to status
	 */
	int IsPerfect(uint genotype, const BitSetType &mask);
	
	 /**
	  * @brief Returns the unique ID for this particular snp	
	  */
	uint GetID();

	/**
	 * @brief Concatinates the contents of other to the local one. 
	 * This is intended for cases where penotypes are broken into separate files
	 */
	bool Merge(const SnpAligned *other);

	/**
	 * @brief Right now, these here to allow them to be sorted in STL containers. 
	 * They don't  perform any sort of data oriented evaluations
	 */
	bool operator==(const SnpAligned& other);
	bool operator<(const SnpAligned& other);

	int ReduceInstanceCount();			///<Used to help with the lifetime within one or more pools
	int IncrementInstanceCount();		///<Used to help with the lifetime within one or more pools

	/**
	 * @brief Label a snp with information associated the order it was loaded from disk
	 * @param labelCount The number of labels passed
	 * @param labels The array of values associated with the current snp
	 * @note A Snp Has multiple labels if it is the product of a crossing (such as a 2 snp model)
	 */
	void SetLabel(uint labelCount, Label labels[]);

	void SetBinaryGenotypes(unsigned char encGenotypes);

	/**
	 * @brief Returns the text representation of the snp's label
	 */
	const char *GetLabel();

	string GetLineNumber();

	/**
	 * @brief Returns  a text descriptor for the local snp
	 */
	string GetTxtDescriptor();

	/**
	 * @brief Returns a binary descriptor for the local snp	 (for now a 32 integer)
	 * @note This might not apply to apply snps that could possibly be derived- This still needs
	 * serious consideration
	 */
	uint GetBinDescriptor();

	/**
	 * Performs a few tests to validate the contents of the snp
	 */
	bool Validate();


	string toString();


	/**
	 * @brief Set the genotype for a given individual to value
	 * @param indivId The individual you wish to set. This is a zero based index
	 * @param value The genotype to be set: 1-Homozygous 2-HeteroZygous 3-Homozygous 0-unknown
	 */
	bool SetValue(int indivId, uint value);

	void SetLabel(const Label &label);

	void ResetEvaluations();
	uint GetEvaluationCount();
	void AppendEvaluation(float score);
	float GetEvaluation(uint foldID);

	/**
	 * @brief generate the base gentoype labels.
	 * This currently only handles AA, Aa and aa as well as unknown. It should be called
	 * after the label has been SNP's setup. 
	 */
	void SetupGenotypeLabels();

	/**
	 * @brief Returns the actual genotype label for reporting purposes. 
	 */
	string GetGenotypeLabel(uint gt);
	void SetGenotypeLabel(uint gt, string& label);

	uint GetBitsetCount() { return bitsetCount;}

	/**
	 * @brief Resets (empties) the bitsets
	 */ 
	void Reset();

	BitSetType Verify();

	/**
	 * @brief Returns the high risk cells associated with the model
	 */
	BitSetType GetHrCells();
	void SetHrCells(BitSetType &hrCells);
protected:
	SnpAligned &operator=(const SnpAligned& other);
	/**
	 * Creation is allowed only for the pool
	 */
    SnpAligned(int genotypeCount, uint id);

	SnpAligned(const SnpAligned& other);

	/**
	 * @brief Initialize the bit vectors to the correct size
	 */
	void InitBitsets(uint individuals);

	/**
	 * @brief Allow a snp to be reinitialized. This simply resizes the bitvector array and clears the bits 
	 */
	void Init(uint genotypeCount);
	void Init(uint genotypeCount, uint individuals);
	/**
	 * @brief Destruction is allowed only for the pool
	 */
    ~SnpAligned();

	/**
     * @brief Contains a bitset for each of the genotypes (possible variables)
	 * For each genotype/discrete variable, there will be a bitset whose values represent individuals who
	 * have this genotype at a given SNP. 
	 */
	BitSetType *bitsets;
	uint bitsetCount;		///<This is just the size of the bitset array. 
	string *bitsetLabels;

	/**
	 * @brief Used to determine how many bitsets are required. For genetic data, maxGenotypes will be 3. However, 
	 * if there is a need for expressing looking for linkage to environmental conditions, there is no hard coded limit.
	 */
	uint genotypeCount;	
	
	/**
	 * @brief Count of the individuals contained within the samples used for this snp
	 */
	uint individualCount;

	/**
	 * @brief This is the human readable label which would appear as you might expect for a report. 
	 * If there is no meaningful string representation, it is just the index as a string. Where the string 
	 * and the binary form differ is the portrayal of crossed SNPs. The string version will have an x 
	 * between each component of the label
	 */
	string label;

	/**
	 * @brief In general, this is a numerical label which can be used to associate a SNP within a set of results.
	 * If the SNP happens to be the result of crossing two or more SNPs, the result is an array of unsigned integers
	 * where labeln[0] is the leftmost portion, and labeln[n] is the rightmost
	 */
	
	Label *binaryLabel;

 	/**
	 * @brief The number of binary labels that can be found
	 */
	uint labelCount;

	/**
	 * @brief Simple integer to allow the pool to identify what is what
	 * @note This isn't a good variable to use for reporting, since multiple SNPs could be encountered with the same
	 * id. Use the label for reporting instead
	 */
	uint id;

	vector<float> mdEvals;

	/**
	 * @brief Used to help keep track of when an item is ready to be moved into the available queue
	 */
	int _instanceCount;

	/**
	 * @brief Converter used to help interpret the meaning of the genetic pairs as they come in
	 * @note This is used ONLY if 
	 */
	GenoBPLookup genLookup;

	/**
	 * @brief Used for logging snps. 
	 * Set this to INCLUDED or a reason for it's exclusion
	 */
	string inclusionStatus;
	
	bool hrCellsSet;
	BitSetType hrCells;
};

inline
BitSetType SnpAligned::GetHrCells() {
	assert(hrCellsSet);
	return hrCells;
}

inline
void SnpAligned::SetHrCells(BitSetType &cells) {
	hrCellsSet = true;
	hrCells = cells;
}

inline
void SnpAligned::ResetEvaluations() {
	mdEvals.clear();
}

inline
uint SnpAligned::GetEvaluationCount() {
	return mdEvals.size();
}	

inline
void SnpAligned::AppendEvaluation(float score) {
	mdEvals.push_back(score);
}

inline
float SnpAligned::GetEvaluation(uint foldID) {
	assert(foldID < mdEvals.size());
	return mdEvals[foldID];
}


/**
 * @brief Quick self test to id any missing data that isn't properly recorded
 */
inline
BitSetType SnpAligned::Verify() {
	BitSetType dataReallyPresent;
	dataReallyPresent = bitsets[1];
	
	for (uint i=2; i<genotypeCount; i++) 
		dataReallyPresent |= bitsets[i];

	BitSetType dataReallyMissing = bitsets[0] - dataReallyPresent;
	return dataReallyMissing;

}

inline
bool SnpAligned::Validate() {
	inclusionStatus="INCLUDED";
	int representedGenotypes = 0;

	for (uint i=1; i<genotypeCount; i++) {
		if (CountIndividuals( i) > 0)
			representedGenotypes++;
	}
	
	if (representedGenotypes==0)
		inclusionStatus="EXCLUDED - NO SNP DATA";
	else if (representedGenotypes==1)
		inclusionStatus="EXCLUDED - CONSTANT GENOTYPE";
	return representedGenotypes>1;
}

inline
string SnpAligned::GetLineNumber() {
	char lineNumber[2048];
	sprintf(lineNumber, "%d", binaryLabel[0].lineNumber);
	
	for (uint i=1; i<labelCount; i++)
		sprintf(lineNumber, "%sx%d", lineNumber, binaryLabel[i].lineNumber);

	return lineNumber;
}


inline
const char *SnpAligned::GetInclusionStatus() {
	return inclusionStatus.c_str();
}

inline
uint SnpAligned::GetLabelCount() {
	return labelCount;
}

inline
uint SnpAligned::GetLabel(uint labelidx) {
	return binaryLabel[labelidx].label;
}

inline
uint SnpAligned::GetLineNumber(uint labelidx) {
	return binaryLabel[labelidx].lineNumber;
}

inline
void SnpAligned::SetBinaryGenotypes(unsigned char encGenotypes) {
	genLookup.ParseBinaryGenotypes(encGenotypes);
}

inline
string SnpAligned::GetTxtDescriptor() {
	return label+genLookup.GetHeterozygote();
}

inline
uint SnpAligned::GetBinDescriptor() {
	//This just looks bad! 8 bits is too small to hold 6 + 3 + 3 bits (chromosome + allele1 + allele2)
	unsigned char desc=binaryLabel[0].label&63;
	desc<<=3;
	desc=desc|genLookup.GetEncodedGenotypes();
	assert(0);
}


inline
void SnpAligned::SetLabel(const Label &label) {
	stringstream lbl;

	if (binaryLabel)
		delete[] binaryLabel;
	binaryLabel=new Label[1];

	binaryLabel[0]=label;
	labelCount=1;
	lbl<<label.label;
	this->label=lbl.str();
}

inline
void SnpAligned::SetLabel(uint labelCount, Label labels[] ) {
	stringstream lbl;
	string sep="x";
	
	if (binaryLabel)
		delete[] binaryLabel;
	
	binaryLabel = new Label[labelCount];
	lbl<<labels[0].label;
	binaryLabel[0]=labels[0];
	for (uint i=1; i<labelCount; i++) {
		binaryLabel[i]=labels[i];
		lbl<<sep<<labels[i].label;
	}
	string s = lbl.str();
	label=s;
	this->labelCount=labelCount;
}
inline
const char *SnpAligned::GetLabel() {
	return label.c_str();
}

inline
int SnpAligned::ReduceInstanceCount() {
	return --_instanceCount;
}
inline
int SnpAligned::IncrementInstanceCount() {
	return ++_instanceCount;
}

inline
SnpAligned &SnpAligned::operator=(const SnpAligned& o) {

	//Having issues with the bitset copy- But, the pool is better for cleaning up copies, so, this isn't required
	assert(0);
	genotypeCount=o.genotypeCount;
	individualCount=o.individualCount;
	label=o.label;
	genLookup=o.genLookup;
		
	for (uint i=0; i<genotypeCount; i++) 
		bitsets[i]=o.bitsets[i];

	return *this;
}

inline
uint SnpAligned::GetID() {	
	return id;
}

inline
bool SnpAligned::operator==(const SnpAligned& other) {
	return id == other.id;
}

inline
bool SnpAligned::operator<(const SnpAligned& other) {
	return id < other.id;
}




inline
void SnpAligned::SetIndividualCount(int individualCount) {
	this->individualCount=individualCount;
	//Clear and resize the bitset
	InitBitsets(individualCount);	

}


inline
uint SnpAligned::CountIndividuals(uint genotypeID) {
	if (genotypeID >= genotypeCount)
		cout<<"Trying to squeeze "<<genotypeID<<" genotypes out of a SNP with only "<<genotypeCount<<" genotypes\n";
	
	assert(genotypeID >=0 && genotypeID < (uint)genotypeCount);

	return bitsets[genotypeID].count();
}


inline
uint SnpAligned::CountGenotypes() {
	return genotypeCount;
}

inline
string SnpAligned::GetGenotypeAlleles(uint genotypeID) {
	assert(genotypeID >=0 && genotypeID<genotypeCount);
	
	//For now, this is a stupid approach 
	if (genotypeID == 0)
		return genLookup.GetHomozygote1();
	else if (genotypeID == 1)
		return genLookup.GetHeterozygote();
	else if (genotypeID == 2)
	 	return genLookup.GetHomozygote2();
	return "MM";
}

inline
BitSetType &SnpAligned::GetGenotype(uint genotypeID) {
	assert(genotypeID >=0 && genotypeID <genotypeCount);

	return bitsets[genotypeID];
}

inline
SnpAligned::SnpAligned(int genotypeCount, uint id) : genotypeCount(genotypeCount), individualCount(0), binaryLabel(NULL), id(id), _instanceCount(0) {
	bitsetCount = genotypeCount;
	bitsets = new BitSetType[genotypeCount];
	bitsetLabels = new string[genotypeCount];
	genLookup.SetDelimeter('\t');
}


inline
void SnpAligned::Init(uint genotypeCount, uint individuals) {
	hrCellsSet = false;
	genLookup.Reset();
	if (genotypeCount > bitsetCount)	{
		this->genotypeCount = bitsetCount = genotypeCount;
		
		if (bitsets)
			delete[] bitsets;
		bitsets = new BitSetType[genotypeCount];
		if (bitsetLabels)
			delete[] bitsetLabels;
		bitsetLabels = new string[genotypeCount];
	}
	this->genotypeCount = genotypeCount;
	if (individuals != bitsets[0].size()) {
		bitsets[0].resize(individuals, true);
		for (uint i=1; i<genotypeCount; i++) 
			bitsets[i].resize(individuals, false);
	}
}



inline
void SnpAligned::SetupGenotypeLabels() {
	static string labelSets[] = {"?/?", "1/1", "1/2", "2/2" }; 
	assert(genotypeCount <= 4);

	for (uint i=0; i<genotypeCount; i++) {
		bitsetLabels[i] = labelSets[i];
	}
}

inline
void SnpAligned::SetGenotypeLabel(uint gt, string& label){
	assert(genotypeCount);
	if (gt > genotypeCount) {
		cout<<"Trying to resize lables beyond size\n";
		abort();
	}
	bitsetLabels[gt] = label;
}

inline
string SnpAligned::GetGenotypeLabel(uint gt) {
	assert(genotypeCount);
	if (gt > genotypeCount) {
		cout<<"Trying to resize lables beyond size\n";
		abort();
	}
	return bitsetLabels[gt];
}



inline
void SnpAligned::Init(uint genotypeCount) {
	hrCellsSet = false;
	genLookup.Reset();
	if (genotypeCount < bitsetCount)	{
		bitsetCount = genotypeCount;
		if (bitsets)
			delete[] bitsets;
		bitsets = new BitSetType[genotypeCount];
		if (bitsetLabels)
			delete[] bitsetLabels;
		bitsetLabels = new string[genotypeCount]; 
		this->genotypeCount = genotypeCount;
		
	}
	this->genotypeCount = genotypeCount;
	
	/**
	 * Reset the bitvectors to contain only 0s
	 */
	bitsets[0].set();
	for (uint i = 1; i<genotypeCount; i++)
		bitsets[i].reset();
}


}

#endif
