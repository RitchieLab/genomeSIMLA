/* 
 * File:   feature.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 11:52 AM
 */

#include "feature.h"
#include "magicdb.h"
#include <sstream>
#include <iomanip>
namespace Paris {

bool Feature::IgnorePValueOfZero = false;

Feature::Feature(uint id, const char* chr, uint beg, uint end) : GenomicRegion(id, chr, beg, end), binIndex((uint)-1), sigCount() { }
Feature::Feature() : binIndex((uint)-1), sigCount(0) { }
Feature::~Feature() { }

bool Feature::AddValue(uint snpIndex, const char *chrom, uint position, float pvalue) {
	if (pvalue > 0 || !IgnorePValueOfZero) {
		if (chrom == _chromosome && position >= _begin && position <= _end) {
			if (ParisResults::resultsDB)
				ParisResults::db.sociDB<<"INSERT INTO feature_snps VALUES (:featureID, :rs, :pos, :pvalue)", use(id), use(snpIndex), use(position), use((double)pvalue);
			
			//assert(pscores.find(position) == pscores.end());
			pscores[position] = pvalue;
			//pvalues.insert(std::pair<float, uint>(pvalue, position));
			if (pvalue > 0 && pvalue <= pvThreshold)
				sigCount++;
			return true;
		}
	}
	return false;
}
void Feature::WriteBinReport(std::ostream& os) {
	std::map<uint, float>::iterator itr = pscores.begin();
	std::map<uint, float>::iterator end = pscores.end();

	float minScore = 2.0;

	while (itr != end) {
		if (itr->second < minScore)
			minScore = itr->second;
		itr++;
	}
	if (minScore > 0.0 and minScore < 2.0)
		os<<binIndex<<"\t"<<minScore<<"\n";
}
void Feature::DetailedReport(std::map<uint, uint>& snps, const char *prefix, std::ostream& os, uint& totalSig, uint &totalNSig) {
	std::stringstream details;
	details<<prefix<<","<<_begin<<","<<_end;

	std::map<uint, float>::iterator itr = pscores.begin();
	std::map<uint, float>::iterator end = pscores.end();
	uint sCount = 0, nsCount = 0;
	while (itr != end) {
		if (itr->second < pvThreshold && itr->second > 0.00000){
			sCount++;
			os<<BinIndex()<<","<<pscores.size()<<","<<
				 
				 
				 details.str()<<",1,rs"<<snps[itr->first]<<","<<itr->first<<","<<itr->second<<"\n";
		}
		else {
			os<<BinIndex()<<","<<pscores.size()<<","<<
				 details.str()<<",0,rs"<<snps[itr->first]<<","<<itr->first<<","<<itr->second<<"\n";
			nsCount++;
		}
		itr++;
	}
	totalSig+=sCount;
	totalNSig+=nsCount;
}

std::map<uint, float> Feature::GetPValues() {
	return pscores;
}
int Feature::CountSignificantMembers() const {
	return sigCount;
}

std::set<uint> Feature::GeneIDs() {
	return geneIDs;
}

uint Feature::BinIndex() {
	return binIndex;
}

void Feature::BinIndex(uint idx) {
	binIndex = idx;
}

uint Feature::FeatureSize() const {
	return pscores.size();
	//return pvalues.size();
}



#ifdef TEST_APP

TEST(FeatureTest, AddValues) {
	Feature::IgnorePValueOfZero = true;
	Feature f(1, "1", 100, 201);
	EXPECT_EQ(true, f.AddValue(100, "1", 101, 0.1));
	EXPECT_EQ(false, f.AddValue(101,"1", 99, 0.1));
	EXPECT_EQ(false, f.AddValue(102,"2", 101, 0.1));
	f.AddValue(103, "1", 110, 0.05);
	EXPECT_EQ(false, f.AddValue(104, "1", 120, 0.0));
	EXPECT_EQ(2, f.FeatureSize());

}


TEST(FeatureTest, Significance) {
	Feature f(1, "1", 100, 201);
	f.AddValue(100, "1", 101, 0.1);
	f.AddValue(101, "1", 150, 0.01);
	f.AddValue(102, "1", 155, 0.0);
	f.AddValue(103, "1", 160, 0.05);

	EXPECT_EQ(2, f.CountSignificantMembers());
}

TEST(FeatureTest, PValues) {
	Feature::IgnorePValueOfZero = true;
	Feature f(1, "1", 100, 201);
	f.AddValue(100, "1", 101, 0.1);
	f.AddValue(101, "1", 150, 0.01);
	f.AddValue(102, "1", 155, 0.0);
	f.AddValue(103, "1", 160, 0.05);
	f.AddValue(104, "1", 165, 0.05);

	std::map<uint, float> pvorig = f.GetPValues();
	std::multimap<float, uint> pvalues;
	std::map<uint, float>::iterator itr = pvorig.begin();
	std::map<uint, float>::iterator end = pvorig.end();

	while (itr != end) {
		pvalues.insert(std::pair<float, uint>(itr->second, itr->first));
		itr++;
	}
	EXPECT_EQ(2, pvalues.count(0.05));
	EXPECT_EQ(1, pvalues.count(0.1));
	EXPECT_EQ(0, pvalues.count(0.0));
	EXPECT_EQ(0, pvalues.count(1.0));
}

#endif //TEST_APP
}
