/* 
 * File:   gene.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 5:01 PM
 */

#include "gene.h"
#include <iostream>
#include "magicdb.h"
#include <sstream>
using namespace std;

namespace Paris {

bool Gene::AllowRedundantFeatures = false;

Gene::Gene() {}
Gene::Gene(const char *ensID, uint id, const char *chr, uint beg, uint end) : GenomicRegion(id, chr, beg, end), ensemblID(ensID), alias(ensID) { }
Gene::Gene(const Gene& other) : ensemblID(other.ensemblID) {}
Gene::~Gene() {}

std::string Gene::EnsemblID(){
	return ensemblID;
}

std::string Gene::Alias() {
	return alias;
}

std::multimap<uint, uint> Gene::GetGroups() {
	return groups;
}

void Gene::SetAlias(const char *alias) {
	this->alias = alias;
}

void Gene::AddFeature(Feature* feature) {
	
	//Four possible ways of being a part of the gene
	// |--{---|    }
	//    {  |--|  }
	//    {   |----}----|
	// |--{--------}----|
	if (feature->_begin <= _end && feature->_end >= _begin) {
		assert(feature->FeatureSize() != (uint)-1);
		if (ParisResults::resultsDB)
			ParisResults::db.sociDB<<"INSERT INTO gene_to_feature VALUES (:geneID, :featureID)", use(id), use(feature->id);
		features.insert(feature);
	}
}

void Gene::DetailedReport(std::map<uint, uint>& snps, const char *prefix, std::ostream& os) {
	std::stringstream details;
	details<<prefix<<","<<id<<","<<_chromosome<<","<<_begin<<","<<_end;
	uint sigPValues = 0, nsigPValues=0;
	std::set<Feature*>::const_iterator itr = features.begin();
	std::set<Feature*>::const_iterator end = features.end();
	while (itr != end) {
		(*itr++)->DetailedReport(snps, details.str().c_str(), os, sigPValues, nsigPValues);
	}
}

std::set<Feature*> Gene::GetFeatures() {
	return features;
}

void Gene::GetFeatureMap(std::multiset<uint>& bins, std::set<uint>& featureIDs) {
	std::set<Feature*>::const_iterator itr = features.begin();
	std::set<Feature*>::const_iterator end = features.end();

	if (features.size() > 0) {
		while (itr != end) {
			Feature *f = (*itr);
			assert(f->FeatureSize() != (uint)-1);
			if (AllowRedundantFeatures || featureIDs.find(f->id) == featureIDs.end()) {
				if (f->FeatureSize() > 0 && f->BinIndex() != (uint)-1) {
					featureIDs.insert(f->id);
					uint  binID = f->BinIndex();
					bins.insert(binID);
				}
			}
			itr++;
		}
	}
}

void Gene::ComprehensiveReport(std::ostream& os) {
	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end) {
		Feature *feature = *itr++;

		os<<"\t"<<id<<"\t"<<EnsemblID()<<"\t"<<_chromosome<<"\t"<<_end-_begin<<"\t"<<"\t"<<feature->_end-_begin<<"\t"<<"\t"<<feature->FeatureSize()<<"\t"<<feature->CountSignificantMembers()<<"\n";

	}
}

void Gene::CollectSNPs(std::set<uint>& snps) {
	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end) {
		Feature *feature = *itr++;
		std::map<uint, float> pvalues = feature->GetPValues();
		std::map<uint, float>::iterator pitr = pvalues.begin();
		std::map<uint, float>::iterator pend = pvalues.end();

		while (pitr != pend) 
			snps.insert(pitr++->first);
	}
}
void Gene::CollectPValues(vector<float>& pvalues) const {
	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	std::set<uint> snps;							///< We'll keep SNPs we've added in here, to avoid repeating them when they appear in multiple features
	while (itr != end ) {
		Feature *feature = *itr++;

		std::map<uint, float> pv = feature->GetPValues();

		std::map<uint, float>::iterator pitr = pv.begin();
		std::map<uint, float>::iterator pend = pv.end();

		while (pitr != pend) {
			if (snps.find(pitr->first) == snps.end()) {
				snps.insert(pitr->first);
				pvalues.push_back(pitr->second);
			}
			pitr++;
		}
	}
}

void Gene::ListFeatures(ostream& os) {
	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	os<<"Gene: "<<id<<"\n";
	int i=0;
	while (itr != end) {
		os<<(i++) + 1<<"\t"<<(*itr)->id<<"\t"<<(*itr)->_begin<<"\t"<<(*itr)->_end<<"\t"<<(*itr)->GetPValues().size()<<"\n";
		itr++;
	}
}

int Gene::CountSignificantMembers(std::set<Feature*>& usedFeatures) const {
	int count = 0;

	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end ) {
		Feature *feature = *itr++;

		if (AllowRedundantFeatures || usedFeatures.find(feature) == usedFeatures.end()) {
			if (feature->CountSignificantMembers() > 0) {
				count++;
			}
			usedFeatures.insert(feature);
		}
	}
	return count;
}

int Gene::CountSignificantMembers(std::set<Feature*>& usedFeatures, uint& simpleFeature, uint &complexFeature, uint& sigSimple, uint &sigComplex) const {
	int count = 0;

	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end ) {
		Feature *feature = *itr++;

		if (AllowRedundantFeatures || usedFeatures.find(feature) == usedFeatures.end()) {
			if (feature->FeatureSize() == 1)
				simpleFeature++;
			else
				complexFeature++;
			if (feature->CountSignificantMembers() > 0) {
				if (feature->FeatureSize() == 1)
					sigSimple++;
				else
					sigComplex++;
				count++;
			}
			usedFeatures.insert(feature);
		}
	}
	return count;
}
int Gene::CountComplexFeatures() {
	int count = 0;

	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end ) {
		if ((*itr++)->FeatureSize() > 1)
			 count++;
	}
	return count;

}

int Gene::CountSimpleFeatures() {
	int count = 0;

	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end ) {
		if ((*itr++)->FeatureSize() == 1)
			 count++;
	}
	return count;
}


int Gene::CountEmptyFeatures() {
	int count = 0;

	std::set<Feature*>::iterator itr = features.begin();
	std::set<Feature*>::iterator end = features.end();

	while (itr != end ) {
		if ((*itr++)->FeatureSize() == 0)
			 count++;
	}
	return count;
}

uint Gene::FeatureCount() const {

	return features.size();
}



void Gene::Summarize(std::ostream& os, bool html) {
	uint simple=0, simpleSig=0, complex=0, complexSig=0;
	std::set<Feature*> usedFeatures;
	CountSignificantMembers(usedFeatures, simple, complex, simpleSig, complexSig);
	if (html)
		os<<"<TR><TD>"<<id<<"</TD><TD><A HREF='http://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g='"
			<<ensemblID<<"'>"<<ensemblID<<"</A></TD><TD>TBD</TD><TD>"<<complex<<"</TD><TD>"<<complexSig<<"</TD><TD>"
			<<simple<<"</TD><TD>"<<simpleSig<<"</TD></TR>";
	else
		os<<id<<","<<ensemblID<<",TBD,"<<complex<<","<<complexSig<<","<<simple<<","<<simpleSig<<"\n";
}
void Gene::AddGroup(uint groupType, uint groupID) {
	groups.insert(pair<uint, uint>(groupType, groupID));
	if (ParisResults::resultsDB)
		ParisResults::db.sociDB<<"INSERT INTO pathway_to_gene VALUES (:groupID, :geneID)", use(groupID), use(id);

}


#ifdef TEST_APP

TEST(GeneTest, AddFeature) {
	Gene gene("ENSG001", 1, "1", 100, 1000);
	Feature f1(1, "1", 75, 125);
	f1.AddValue(1, "1", 80, 0.01);
	f1.AddValue(2, "1", 85, 0.40);

	Feature f2(2, "1", 600, 750);
	f2.AddValue(3, "1",630, 0.1);
	f2.AddValue(7, "1",700, 0.01);
	f2.AddValue(8, "1",725, 0.99);

	Feature f3(3, "1", 700, 1225);
	f3.AddValue(4, "1", 751, 0.01);

	Feature f4(4, "1", 75, 1225);
	f4.AddValue(5, "1", 80, 0.01);

	Feature f5(5, "1", 130, 130);
	f5.AddValue(5, "1", 130, 0.001);
	
	Feature f6(6, "1", 550, 550);
	f6.AddValue(6, "1", 550, 0.90);

	gene.AddFeature(&f1);
	gene.AddFeature(&f2);
	gene.AddFeature(&f3);
	gene.AddFeature(&f4);
	gene.AddFeature(&f5);
	gene.AddFeature(&f6);
	EXPECT_EQ(6, gene.FeatureCount());
	EXPECT_EQ(5, gene.CountSignificantMembers());
	uint simpleFeature=0, sigSimple=0, complexFeature=0, sigComplex=0;
	gene.CountSignificantMembers(simpleFeature, complexFeature, sigSimple, sigComplex);
	EXPECT_EQ(4, simpleFeature);
	EXPECT_EQ(2, complexFeature);
	EXPECT_EQ(3, sigSimple);
	EXPECT_EQ(2, sigComplex);

}

TEST(GeneTest, TestGetFeatureMap) {
	Gene gene("ENSG001", 1, "1", 100, 1000);
	Feature f1(1, "1", 75, 125);
	f1.AddValue(1, "1", 80, 0.01);
	f1.AddValue(2, "1", 80, 0.40);
	f1.BinIndex(1);

	Feature f2(2, "1", 600, 750);
	f2.AddValue(3, "1",630, 0.1);
	f2.AddValue(7, "1",700, 0.01);
	f2.AddValue(8, "1",725, 0.99);
	f2.BinIndex(1);

	Feature f3(3, "1", 700, 1225);
	f3.AddValue(4, "1", 751, 0.01);
	f3.BinIndex(2);

	Feature f4(4, "1", 75, 1225);
	f4.AddValue(5, "1", 80, 0.01);
	f4.BinIndex(2);

	Feature f5(5, "1", 130, 130);
	f5.AddValue(5, "1", 130, 0.001);
	f5.BinIndex(2);

	Feature f6(6, "1", 550, 550);
	f6.AddValue(6, "1", 550, 0.90);
	f6.BinIndex(2);

	gene.AddFeature(&f1);
	gene.AddFeature(&f2);
	gene.AddFeature(&f3);
	gene.AddFeature(&f4);
	gene.AddFeature(&f5);
	gene.AddFeature(&f6);


	std::multiset<uint> bins;
	std::set<uint> featureIDs;
	gene.GetFeatureMap(bins, featureIDs);
	EXPECT_EQ(2, bins.count(1));
	EXPECT_EQ(4, bins.count(2));
	EXPECT_EQ(6, featureIDs.size());
	
}

#endif //TEST_APP

}
