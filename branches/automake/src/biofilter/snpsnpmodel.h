//
// C++ Implementation: snpsnpmodel.h
//
// Description: Defines functionality associated with Snp/Snp models
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __SNP_SNP_MODEL_H
#define __SNP_SNP_MODEL_H
#include "biomodel.h"
#include <set>
#include <iomanip>
#include <assert.h>
#ifndef uint 
typedef unsigned int uint;
#endif
namespace Biofilter {



class SnpSnpModel : public BioModel {
public:
	SnpSnpModel(uint snp1, uint snp2, float implicationIndex);
	virtual ~SnpSnpModel();

	/* file buffer requirements */
	bool operator<(const SnpSnpModel& other) const;
	bool EvalLT(const SnpSnpModel* other) const;
	bool operator==(const SnpSnpModel& other) const;
	SnpSnpModel& operator=(const SnpSnpModel& other);
	void WriteBinary(std::ostream& file) const;
	void Write(std::ostream& os, bool useBinary = true) const;
	bool LoadBinary(std::istream& file);
	//void MergeGroups(const SnpSnpModel& other);

	std::set<uint> snps;
	//float implicationIndex;		///< This is really a genegene property....
};

class SnpSnpModel;
struct LtSnpModelPointer {
	bool operator()(const SnpSnpModel *model1, const SnpSnpModel *model2) const {
		return *model1 < *model2;
		//return model1->EvalLT(model2);
	}
};
typedef std::set<SnpSnpModel*, LtSnpModelPointer> SnpModelCollection;

inline
SnpSnpModel::SnpSnpModel(uint snp1, uint snp2, float implicationIndex) : BioModel(implicationIndex) {
	snps.insert(snp1);
	snps.insert(snp2);
	//assert(snp1<10000000);		//This is a test to see if the SNPs are coming in correctly or being broken someplace else
}

inline
SnpSnpModel::~SnpSnpModel() { }

inline
bool SnpSnpModel::operator<(const SnpSnpModel& other) const {
	return snps < other.snps;
}

inline
bool SnpSnpModel::EvalLT(const SnpSnpModel* other) const {
	return snps < other->snps;
}

inline
bool SnpSnpModel::operator==(const SnpSnpModel& other) const {
	return snps == other.snps;
}

inline
SnpSnpModel& SnpSnpModel::operator=(const SnpSnpModel& other) {
	implicationIndex = other.implicationIndex;
	snps = other.snps;
	return *this;
}


inline
void SnpSnpModel::Write(std::ostream& os, bool useBinary) const {
	if (useBinary) {
		WriteBinary(os);
		return;
	}
	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();
	//uint count = 0;
	while (itr != end) {
		os<<std::setw(11)<<std::right<<*itr++;
	}
	os<<std::setw(7)<<std::right<<std::setprecision(1)<<ImplicationIndex()<<"\n";
}

inline
void SnpSnpModel::WriteBinary(std::ostream& file) const {
	//uint modelSize = snps.size();
	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();
	
	while (itr != end) {
		uint locus = *itr++;
		file.write((char*)&locus, sizeof(uint));
	}
	file.write((char*)&implicationIndex, sizeof(float));           //# implication index
}

inline
bool SnpSnpModel::LoadBinary(std::istream& file) {
	uint modelSize = 2;
	//file.read((char*)&modelSize, sizeof(uint));

	for (size_t i=0; i<modelSize; i++) {
		uint locus = 0;
		file.read((char*)&locus, sizeof(uint));
		snps.insert(locus);
	}
	file.read((char*)&implicationIndex, sizeof(float));
	return modelSize > 0;
}

/*
inline
void SnpSnpModel::MergeGroups(const SnpSnpModel& other) {

}
*/
}

#endif //__SNP_SNP_MODEL_H
