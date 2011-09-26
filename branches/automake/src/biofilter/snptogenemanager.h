//
// C++ Interface: snptogenemanager
//
// Description: Simple container object to generate lists of genes
//              according to the rsID and chromosome
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>,
// (C) Marylyn Ritchie 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SNP_TO_GENE_MANAGER
#define SNP_TO_GENE_MANAGER

#include <string>
#include <vector>
#include "utility/strings.h"
#include "kbregion.h"
#include <iostream>
namespace Biofilter {
/**
 * @brief Aid in associating a list of genes to a particular rsID / chromosome combo
 */
class SnpToGeneManager {
public:
	SnpToGeneManager() {}
	~SnpToGeneManager() {}

	/**
	 * @brief Add a gene to the manager object
	 */
	void AddGene(const char *chrom, uint rsID, Knowledge::KbRegion *gene) {
		std::string id = GetID(chrom, rsID);
		snpLists[chrom].insert(rsID);
		geneLists[id].insert(gene);
	}

	/**
	 * @brief Retrieve a list keyed off of the chromosome and rsID
	 */
	std::set<Knowledge::KbRegion*> GetGeneList(const char *chrom, uint rsID) {
		std::set<Knowledge::KbRegion*> geneList;

		if (chrom == NULL) {
			std::map<std::string, std::set<uint> >::iterator chr = snpLists.begin();
			std::map<std::string, std::set<uint> >::iterator chrEnd = snpLists.end();

			while (chr!=chrEnd) {
				std::string id = GetID(chr->first.c_str(), rsID);
				if (geneLists.find(id) != geneLists.end()) {
					std::set<Knowledge::KbRegion*> localList = geneLists[id];
					geneList.insert(localList.begin(), localList.end());
				}
				chr++;
			}
		}
		else {
			std::string id = GetID(chrom, rsID);
			if (geneLists.find(id) != geneLists.end())
				geneList = geneLists[id];
		}
		return geneList;
	}

	std::set<uint> GetSnpLists(const char *chromosome) {
		std::set<uint> snpList;
		if (snpLists.find(chromosome) != snpLists.end())
			snpList = snpLists[chromosome];
		return snpList;
	}

	void PrintReport(std::ostream& os) {
		os<<"Snp To Gene Manager:\n";
		os<<"# Genes: "<<geneLists.size()<<"\n";
		std::map<std::string, std::set<uint> >::iterator itr = snpLists.begin();
		std::map<std::string, std::set<uint> >::iterator end = snpLists.end();

		while (itr != end) {
			os<<itr->first<<"\t"<<itr->second.size()<<"\n";
			itr++;
		}
	}
protected:
	/**
	 * @brief Used to generate the key associated for a given chom/rsid pair
	 */
	std::string GetID(const char *chrom, uint rsID) {
		return std::string(chrom) + "-" + Utility::ToString(rsID);

	}
	std::map<std::string, std::set<Knowledge::KbRegion*> > geneLists;		///< key -> gene
	std::map<std::string, std::set<uint> > snpLists;							///< Chrom->[rsids]
};
}
#endif //SNP_TO_GENE_MANAGER
