/* 
 * File:   chromosome.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 2:35 PM
 */

#include "chromosome.h"
#include <iomanip>
#include <algorithm>
#include "magicdb.h"
#include "parislogs.h"
#include "utility/strings.h"
#include <math.h>

using namespace std;
using namespace soci;

namespace Paris {



Chromosome::Chromosome() : chrID(0), length(0) { }
Chromosome::Chromosome(const char *chrID) : chrID(chrID), length(0) {
	this->chrID.erase(remove_if(this->chrID.begin(), this->chrID.end(), ::isspace), this->chrID.end());
}
Chromosome::Chromosome(const Chromosome& orig) : chrID(orig.chrID), length(orig.length) {
	cout<<"Trying to copy points which shouldn't be copied!\n";
	assert(false);
}
Chromosome::~Chromosome() {
	std::vector<Feature*>::iterator itr = features.begin();
	std::vector<Feature*>::iterator end = features.end();

	while (itr != end)
		delete *itr++;

	std::map<uint, Gene*>::iterator geneItr = genes.begin();
	std::map<uint, Gene*>::iterator endGene = genes.end();

	while (geneItr != endGene)
		delete geneItr++->second;
}



void Chromosome::MergeFeaturesIntoGenes(ostream& os) {
	std::map<uint, Gene*>::iterator geneItr = genes.begin();
	std::map<uint, Gene*>::iterator endGene = genes.end();

	std::set<uint> totalFeatureCounts;
	std::multiset<uint> totalCounts;
	//cerr<<"Chromosome "<<chrID<<"\nChr   Gene      Fefatures\n";
	while (geneItr != endGene) {
		std::set<uint> featureCount;
		std::multiset<uint> counts;
		Gene* gene= geneItr++->second;

		std::set<Feature*> matchedFeatures = GetFeatures(gene->_begin, gene->_end);
		std::set<Feature*>::iterator fItr = matchedFeatures.begin();
		std::set<Feature*>::iterator fEnd = matchedFeatures.end();

		while (fItr != fEnd) {
			Feature *feature = *fItr;
			if (feature->FeatureSize() > 0) {

				//FileLogs::logger.featureDetails<<chrID<<"\t"<<gene->id<<"\t"<<gene->EnsemblID()<<"\t"<<feature->_begin<<"\t"<<feature->_end;
				gene->AddFeature(feature);

				if (FileLogs::WriteFeatureDetails) {
					std::map<uint, float> pvalues = feature->GetPValues();
					std::map<uint, float>::iterator pi = pvalues.begin();
					std::map<uint, float>::iterator pe = pvalues.end();

					while (pi != pe) {
							FileLogs::logger.featureDetails<<chrID<<","<<gene->id<<","<<gene->EnsemblID()<<","<<feature->_begin<<","<<feature->_end<<","<<snps[pi->first]<<","<<pi->first<<","<<pi->second<<"\n";
						//FileLogs::logger.featureDetails<<"\t[ "<<pi->first<<" "<<pi->second<<" ]";
						pi++;
					}
				}
				//if (feature->_begin == feature->_end)
				//	FileLogs::logger.featureDetails<<"\t Singular";
				//FileLogs::logger.featureDetails<<"\n";
			}
/*			else {
				FileLogs::logger.emptyFeatureDetails<<chrID<<"\t"<<gene->id<<"\t"<<gene->EnsemblID()<<"\t"<<feature->_begin<<"\t"<<feature->_end;
				if (feature->_begin == feature->_end)
					FileLogs::logger.emptyFeatureDetails<<"\t Singular";
				FileLogs::logger.emptyFeatureDetails<<"\n";
			}
*/			fItr++;
		}
	}

	//cerr<<"\tChromosome "<<chrID<<" \tGenes: "<<genes.size()<<"\tFeatures: "<<features.size()<<"\n";
}

Gene *Chromosome::GetGene(uint geneID) {
	if (genes.find(geneID) != genes.end())
		return genes[geneID];
	else
		return NULL;
}

void Chromosome::AddGene(const char *ensID, uint id, uint start, uint end, map<uint, string>& aliasLookup) {
	Gene *newGene = new Gene(ensID, id, chrID.c_str(), start, end);

	if (aliasLookup.find(id) != aliasLookup.end()) 
		newGene->SetAlias(aliasLookup[id].c_str());
	if (ParisResults::resultsDB)
		ParisResults::db.sociDB<<"INSERT INTO genes VALUES (:geneID, :ensID, :chrID, :start, :end)", use(id), use(std::string(ensID)), use(chrID), use(start), use(end);

	genes[id] = newGene;
}
uint Chromosome::LoadGenes(soci::session& sociDB, uint popID, uint geneExpansion, map<uint, string>& aliasLookup) {
	rowset<row> rs = (sociDB.prepare << "SELECT gene_id, primary_name, chrom, description, start, end FROM regions NATURAL JOIN region_bounds WHERE chrom=:id AND population_id=:popID GROUP BY regions.gene_id", use(chrID), use(popID));
	uint count = 0;
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint geneID = row.get<int>(0);
		string ensemblID  = row.get<string>(1);
		string chrom		= row.get<string>(2);
		//string desc			= row.get<string>(3);
		uint start			= row.get<int>(4);

		//We can't really afford to wrap around into negative space on the uint....very bad
		if (start > geneExpansion)
			start-=geneExpansion;
		else
			start = 0;

		//I think it's OK to let the gene stretch past the end of the chromosome
		uint end				= row.get<int>(5)+geneExpansion;
		AddGene(ensemblID.c_str(), geneID, start, end, aliasLookup);
		count++;
	}
	return count;
}

void Chromosome::AddFeature(uint geneID, uint start, uint stop) {
	assert(start<=stop);
	Feature *feature = new Feature(geneID, chrID.c_str(), start, stop);
	if (ParisResults::resultsDB)
		ParisResults::db.sociDB<<"INSERT INTO features VALUES (:featureID, :chrID, :start, :end)", use(feature->id), use(chrID), use(start), use(stop);
	featureEnd.insert(std::pair<uint, Feature*>(stop, feature));
	featureStart.insert(std::pair<uint, Feature*>(start,feature));

	features.push_back(feature);
}
uint Chromosome::LoadFeatures(soci::session& sociDB, const char *popID, uint &featureID) {
	std::set<uint> availableSNPs;
	std::map<uint, uint>::iterator itr = snps.begin();
	std::map<uint, uint>::iterator end = snps.end();

	while (itr != end) {
		availableSNPs.insert(itr->first);
		itr++;
	}
	rowset<row> rs = (sociDB.prepare << "SELECT DISTINCT start, stop FROM ld_blocks WHERE chromosome=:id AND population_id=:popID", use(chrID), use(string(popID)));
	uint count = 0;

	uint snpCount = 0;
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint start = row.get<int>(0);
		uint stop  = row.get<int>(1);

		if (stop < start) {
			uint t = start;
			start = stop;
			stop = t;
		}

		AddFeature(start, start, stop);
		//Lets remove any SNPs contained within the feature from the pool of available SNPs
		std::map<uint, uint>::iterator itr = snps.lower_bound(start);
		std::map<uint, uint>::iterator end = snps.upper_bound(stop);
		
		while (itr != end) {
			
			snpCount++;
			//Testing against previous results
			availableSNPs.erase(itr->first);
			itr++;
		}
		count++;
	}

	uint featureCount = features.size();
	std::set<uint>::iterator avSNP = availableSNPs.begin();
	std::set<uint>::iterator avSNPend = availableSNPs.end();

	while (avSNP != avSNPend) {
		uint snp = *avSNP++;
		if (featureEnd.find(snp) == featureEnd.end() && featureStart.find(snp) == featureStart.end()) {
			uint pos=snps[snp];
			AddFeature(pos, snp, snp);
			count++;
		}
	}

	
	cerr<<"Chromosome "<<chrID<<" Total SNPs: "<<snps.size()<<" Loaded "<<features.size()<<" Features ("<<featureCount<<")\n";
	return count;
}
void Chromosome::AddValue(uint snpIndex, float pvalue) {
	const char *c = chrID.c_str();
	if (posLookup.find(snpIndex) == posLookup.end())
		return;
	assert(posLookup.find(snpIndex) != posLookup.end());
	uint pos = posLookup[snpIndex];

	uint featureCount = 0;

	std::map<uint, Feature*>::iterator fItr = featureEnd.lower_bound(pos);
	std::map<uint, Feature*>::iterator fEnd = featureEnd.upper_bound(pos);

	if (ParisResults::resultsDB)
		ParisResults::db.sociDB<<"INSERT INTO snps VALUES (:chrID, :rsID, :pos)", use(chrID), use(snpIndex), use(pos);
	
	while (fItr != fEnd) {
		Feature *f = fItr++->second;
		f->AddValue(snpIndex, c, pos, pvalue);
		featureCount++;
	}
	//Just because the end has passed the local value, the beginning might not have
	if (fEnd != featureEnd.end()) {
		Feature *f = fItr++->second;
		if (f->_begin <= pos) {
			f->AddValue(snpIndex, c, pos, pvalue);
			featureCount++;
		}
	}


	fItr = featureStart.lower_bound(pos);
	fEnd = featureStart.upper_bound(pos);
	if (fItr != featureStart.end()) {
		if (fItr == fEnd) {
			if (featureEnd.find(fItr->first) == featureEnd.end()) {
				fItr->second->AddValue(snpIndex, c, pos, pvalue);
				featureCount++;
			}
		}

		while (fItr != fEnd) {
			Feature *f = fItr++->second;
			if (featureEnd.find(fItr->first) == featureEnd.end()) {
				f->AddValue(snpIndex, c, pos, pvalue);
				featureCount++;
			}
		}
	}
	assert(featureCount>0);
}
void Chromosome::VerifyBins() {
	std::vector<Feature*>::iterator itr = features.begin();
	std::vector<Feature*>::iterator end = features.end();

	while (itr != end) {
		Feature *feature = *itr++;
		if (feature->FeatureSize() > 0) {
			if (feature->BinIndex() > 1000) {
				uint sig=0, nsig=0;
				feature->DetailedReport(snps, "prefix", cerr, sig, nsig);
				cerr<<"Feature ID: "<<feature->id<<"\n";
				assert(feature->BinIndex() < 1000);
			}
		}
	}
}
void Chromosome::InitBins(std::set<Feature*, SortByFeatureSize>& bins, std::set<Feature*>& singleFeatureBins, ostream& os) {
	MergeFeaturesIntoGenes(os);


	std::vector<Feature*>::iterator itr = features.begin();
	std::vector<Feature*>::iterator end = features.end();

	uint singles = 0, complexes = 0;
	while (itr != end) {
		Feature *feature = *itr++;
		//Basically, we want to segregate singular features from complex
		if (feature->FeatureSize() > 0) {
			if (feature->FeatureSize() == 1) {
				singles++;
				singleFeatureBins.insert(feature);
			}
			else if (feature->FeatureSize() > 1) {
				complexes++;
				bins.insert(feature);
			}
			else
				cerr<<"Won't you take me to Wonky Town?\n";
		}
	}

}

string Chromosome::ID() {
	return chrID;
}

uint Chromosome::Length() {
	return length;
}

void Chromosome::AddSNP(uint rsID, uint pos) {
	snps[pos] = rsID;
	posLookup[rsID] = pos;

	if (pos > length)
		length = pos;
}

uint Chromosome::SnpCount() {
	return snps.size();
}

uint Chromosome::FeatureCount() {
	return features.size();
}




std::set<Feature*> Chromosome::GetFeatures(uint start, uint stop) {
	std::set<Feature*> matches;
	
	std::map<uint, Feature*>::iterator itr = featureEnd.lower_bound(start);
	std::map<uint, Feature*>::iterator end = featureEnd.upper_bound(stop);
	
	if (itr != featureEnd.end()) {
		while (itr != end) {
			Feature *f = itr++->second;
			matches.insert(f);
		}

		if (end != featureEnd.end()) {
			//We must also get the end
			Feature *f = end->second;
			if (f->_begin<=stop)
				matches.insert(f);
		}
	}
	
	itr = featureStart.lower_bound(start);
	end = featureStart.upper_bound(stop);
	if (itr != featureStart.end()) {
		while (itr != end) {
			Feature *f = itr++->second;
			matches.insert(f);
		}
	}

	return matches;
}


uint Chromosome::GetSNP(uint pos) {
	return snps[pos];
}

std::map<uint, uint> &Chromosome::GetSNPs() {
	return snps;
}

void Chromosome::WriteBinReport(std::ostream& os) {
	std::vector<Feature*>::iterator itr = features.begin();
	std::vector<Feature*>::iterator end = features.end();

	while (itr != end) {
		(*itr)->WriteBinReport(os);
		itr++;
	}
}

std::map<uint, uint> Chromosome::GetSnps(uint start, uint stop) {
	std::map<uint, uint>::iterator itr = snps.lower_bound(start);
	std::map<uint, uint>::iterator end = snps.upper_bound(stop);

	std::map<uint, uint> matches;
	if (itr != snps.end()) {
		while (itr != end) {
			matches[itr->first] = itr->second;
			itr++;
		}
	}
	return matches;
}

#ifdef TEST_APP

class ChromosomeTest : public ::testing::Test {
protected:
	ChromosomeTest() { 
		//ASSERT_EQ(1, Utility::FileExists("SimKB.sqlite"));
		string cnxParam = "dbname=SimKB.sqlite timeout=10";
		sociDB.open(soci::sqlite3, cnxParam.c_str());
		rowset<row> rs = (sociDB.prepare << "SELECT * FROM dataset WHERE chromosome=:id", use(std::string("1")));
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			uint rsid = row.get<int>(1);
			uint pos  = row.get<int>(2);
			string pvs = row.get<string>(3);
			double pvalue = atof(pvs.c_str());
			//double pvalue = row.get<float>(3);
			snps[pos] = rsid;
			pvalues[rsid] = pvalue;
		}
	}
	~ChromosomeTest() {

	}
	virtual void SetUp() {

	}
	virtual void TearDown() {

	}

	void LoadDataset(Chromosome& c) {
		std::map<uint, uint>::iterator itr = snps.begin();
		std::map<uint, uint>::iterator end = snps.end();

		while (itr != end) {
			c.AddSNP(itr->second, itr->first);
			itr++;
		}

	}

	soci::session sociDB;

	std::map<uint, uint> snps;			//pos -> rsid
	std::map<uint, float> pvalues;
};
TEST_F(ChromosomeTest, LoadSNPs) {
	Chromosome c("1");
	LoadDataset(c);

	EXPECT_EQ(snps.size(), c.SnpCount());
	EXPECT_EQ(snps.size(), c.GetSnps(0, c.Length()).size());
	uint count;
	sociDB<<"SELECT COUNT(*) FROM dataset WHERE chromosome='1' AND pos>=:start AND pos<=:stop", use(91000), use(93000), into(count);
	EXPECT_EQ(count, c.GetSnps(91000, 93000).size());

}

TEST_F(ChromosomeTest, LoadFeatures) {
	Chromosome c("1");
	LoadDataset(c);
	uint id;									///< This is used by the load command to assign unique ids
	c.LoadFeatures(sociDB, "CEU", id);
	uint count;
	sociDB<<"SELECT COUNT(*) FROM features WHERE chromosome='1' AND stop>=:begin AND start<=:end", use(91000), use(93000), into(count);
	std::set<Feature*> features	= c.GetFeatures(91000, 93000);
	std::set<Feature*>::iterator fitr = features.begin();
	std::set<Feature*>::iterator fend = features.end();

	while (fitr != fend) {
		Feature *f = *fitr++;
		cerr<<"\t"<<f->_begin<<"\t"<<f->FeatureSize()<<"\n";
	}
	EXPECT_EQ(count, features.size());
	sociDB<<"SELECT COUNT(*) FROM features WHERE chromosome='1' AND stop>=:begin AND start<=:end", use(1500), use(75000), into(count);
	features = c.GetFeatures(1500, 75000);
	while (fitr != fend) {
		Feature *f = *fitr++;
		cerr<<"\t"<<f->_begin<<"\t"<<f->FeatureSize()<<"\n";
	}
	EXPECT_EQ(count, features.size());
}

TEST_F(ChromosomeTest, MergeIntoGenesHandMade) {
	Chromosome c("1");
	                              // 0   1   2   3   4   5
	c.AddSNP(1, 1326000);			//  
	c.AddSNP(2, 1326088);			// X
	c.AddSNP(3, 1326999);			// X   X
	c.AddSNP(4, 1327000);			//     X
	c.AddSNP(5, 1327500);			//     X
	c.AddSNP(6, 1327764);         //     X   X
	c.AddSNP(7, 1327800);         //         X
	c.AddSNP(8, 1328955);         //         X   X
	c.AddSNP(9, 1329000);         //             X
	c.AddSNP(10, 1330000);        //             X
	c.AddSNP(11, 1330065);        //             X
	c.AddSNP(12, 1335500);        // <- Single in gene(5)
	c.AddSNP(13, 1346000);

	c.AddFeature(0, 1326088, 1326999);
	c.AddFeature(1, 1326999, 1327764);
	c.AddFeature(2, 1327764, 1328955);
	c.AddFeature(3, 1328955, 1330065);
	c.AddFeature(4, 1330100, 1335000);
	//Since we aren't using LoadFeatures, we have to add our singles ourselves
	c.AddFeature(5, 1335500, 1335500);


	EXPECT_EQ(2, c.GetFeatures(1326100, 1327700).size());
	EXPECT_EQ(3, c.GetFeatures(1326100, 1327764).size());
	EXPECT_EQ(4, c.GetFeatures(1326100, 1328955).size());

																// 0  1  2  3  4
	c.AddGene("ENSG00000", 0, 1326200, 1326999);	// X  X
	c.AddGene("ENSG00001", 1, 1327000, 1328000);	//    X  X
	c.AddGene("ENSG00002", 2, 1326500, 1329000);	// X  X  X  X
	c.AddGene("ENSG00003", 3, 1325000, 1330000);	// X  X  X  X
	c.AddGene("ENSG00004", 4, 1326999, 1327764);	// X  X  X
	c.AddGene("ENSG00005", 5, 1330065, 1340000);	//             X (plus 1 singular)

											// 0   1   2   3   4
	c.AddValue(1, 0.500);				//
	c.AddValue(2, 0.250);				// X
	c.AddValue(3, 0.040);				// X   X
	c.AddValue(4, 0.010);				//     X
	c.AddValue(5, 0.020);				//     X
	c.AddValue(6, 0.005);				//     X   X
	c.AddValue(7, 0.100);				//         X
	c.AddValue(8, 0.250);				//         X   X
	c.AddValue(9, 0.300);				//             X
	c.AddValue(10, 0.500);				//             X
	c.AddValue(11, 0.500);				//             X
	c.AddValue(12, 0.500);				// <- Single in gene(5)
	c.AddValue(13, 0.500);				//


	std::set<Feature*> features = c.GetFeatures(1327000, 1327100);
	ASSERT_EQ(1, features.size());
	Feature *f = *(features.begin());
	EXPECT_EQ(4, f->FeatureSize());
	EXPECT_EQ(4, f->CountSignificantMembers());

	std::vector<float> pvalues;
	c.MergeFeaturesIntoGenes(cerr);
	Gene* g = c.GetGene(0);
	EXPECT_EQ(2, g->FeatureCount()) << "Feature Count Gene 0\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(5, pvalues.size()) << "PValue Count Gene 0\n";
	pvalues.clear();
	g = c.GetGene(1);
	EXPECT_EQ(2, g->FeatureCount()) << "Feature Count Gene 1\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(6, pvalues.size()) << "PValue Count Gene 1\n";
	pvalues.clear();
	g = c.GetGene(2);
	EXPECT_EQ(4, g->FeatureCount()) << "Feature Count Gene 2\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(10, pvalues.size()) << "PValue Count Gene 2\n";
	pvalues.clear();
	g = c.GetGene(3);
	EXPECT_EQ(4, g->FeatureCount()) << "Feature Count Gene 3\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(10, pvalues.size()) << "PValue Count Gene 3\n";
	pvalues.clear();
	g = c.GetGene(4);
	EXPECT_EQ(3, g->FeatureCount()) << "Feature Count Gene 4\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(7, pvalues.size()) << "PValue Count Gene 4\n";
	pvalues.clear();
	g = c.GetGene(5);
	EXPECT_EQ(2, g->FeatureCount()) << "Feature Count Gene 5\n";
	g->CollectPValues(pvalues);
	EXPECT_EQ(5, pvalues.size()) << "PValue Count Gene 5\n";
}
TEST_F(ChromosomeTest, MergeIntoGenes) {
	Chromosome c("1");
	LoadDataset(c);
	uint id = 0;									///< This is used by the load command to assign unique ids
	c.LoadFeatures(sociDB, "CEU", id);
	c.LoadGenes(sociDB, 0, 0);
	std::map<uint, float>::iterator itr = pvalues.begin();
	std::map<uint, float>::iterator end = pvalues.end();

	while (itr != end) {
		c.AddValue(itr->first, itr->second);
		itr++;
	}

	c.MergeFeaturesIntoGenes(cerr);
	uint count=0, geneCount=0;
	sociDB<<"SELECT COUNT(*) FROM regions", into(geneCount);


	for (uint i=1; i<=geneCount; i++ ) {
		std::vector<float> pvalues;
		sociDB<<"SELECT COUNT(*) FROM (SELECT DISTINCT start FROM features JOIN dataset ON (pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1' AND features.gene_id=:i AND dataset.gene_id=:i)) a", use(i), into(count);
		//sociDB<<"SELECT COUNT(*) FROM (SELECT start,COUNT( rsid) FROM features JOIN dataset ON ( pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1' AND dataset.gene_id=:i AND features.gene_id=:i) GROUP BY start) a", use(i), use(i), into(count);
		//sociDB<<"SELECT COUNT(*) FROM (SELECT start, rsid FROM features JOIN dataset ON (pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1' AND dataset.gene_id=:i AND features.gene_id=:i) GROUP BY start) a", use(i), use(i), into(count);
		//sociDB<<"SELECT COUNT(*) FROM features WHERE gene_id=:geneID", use(i), into(count);
		Gene *gene = c.GetGene(i);
		if (gene) {
			EXPECT_EQ(count, gene->FeatureCount()) << "Gene ("<<i<<") -> Features";
			gene->CollectPValues(pvalues);
			sociDB<<"SELECT COUNT(DISTINCT rsid) FROM features JOIN dataset ON (pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1' AND dataset.gene_id=:i AND features.gene_id=:i)", use(i), into(count);
			if (pvalues.size() != count) {
				cerr<<"SELECT COUNT(DISTINCT rsid) FROM features JOIN dataset ON (pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1' AND dataset.gene_id="<<i<<" AND features.gene_id="<<i<<")\n";
				sort(pvalues.begin(), pvalues.end());
				std::vector<float>::iterator itr = pvalues.begin();
				std::vector<float>::iterator end = pvalues.end();

				while (itr != end) {
					cerr<<"\t"<<*itr++<<"\n";
				}
			}
			EXPECT_EQ(count, pvalues.size()) << "Gene ("<<i<<") -> SNP count";
		}
	}

}

TEST_F(ChromosomeTest, InitBins) {
	Chromosome c("1");
	LoadDataset(c);
	uint id;									///< This is used by the load command to assign unique ids
	c.LoadFeatures(sociDB, "CEU", id);
	c.LoadGenes(sociDB, 0, 0);
	std::map<uint, float>::iterator itr = pvalues.begin();
	std::map<uint, float>::iterator end = pvalues.end();

	while (itr != end) {
		c.AddValue(itr->first, itr->second);
		itr++;
	}

	c.MergeFeaturesIntoGenes(cerr);
	std::set<Feature*, SortByFeatureSize> bins;
	std::set<Feature*> singleFeatureBins;
	c.InitBins(bins, singleFeatureBins, cerr);

	uint geneCount;
	sociDB<<"SELECT COUNT(*) FROM regions WHERE chrom='1'", into(geneCount);

	uint complex, single;
	sociDB<<"SELECT COUNT(DISTINCT start) FROM (SELECT start, COUNT(rsid) AS feature_count FROM features JOIN dataset ON(pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1') GROUP BY start) a WHERE a.feature_count>1", into(complex);
	sociDB<<"SELECT COUNT(DISTINCT start) FROM (SELECT start, COUNT(rsid) AS feature_count FROM features JOIN dataset ON(pos>=start AND pos<=stop AND features.chromosome='1' AND dataset.chromosome='1') GROUP BY start) a WHERE a.feature_count=1", into(single);
	EXPECT_EQ(single, singleFeatureBins.size());
	std::set<Feature*>::iterator bitr = singleFeatureBins.begin();
	std::set<Feature*>::iterator bend = singleFeatureBins.end();


	while (bitr != bend) {
		cerr<<"\t"<<(*bitr)->_begin<<"\t"<<(*bitr)->_end<<"\t"<<(*bitr)->FeatureSize()<<"\n";
		bitr++;
	}

	EXPECT_EQ(complex, bins.size());

	bitr = bins.begin();
	bend = bins.end();

	while (bitr != bend) {
		cerr<<"\t"<<(*bitr)->_begin<<"\t"<<(*bitr)->_end<<"\t"<<(*bitr)->FeatureSize()<<"\n";
		bitr++;
	}

}

#endif //TEST_APP
}
