/* 
 * File:   gene.h
 * Author: torstees
 *
 *
 *
 * Created on January 5, 2010, 5:01 PM
 */

#ifndef _GENE_H
#define _GENE_H
#include "genomicregion.h"
#include "feature.h"
#include <set>
#include <string>
#include <vector>

namespace Paris {

/**
 * Represents a gene which can be a member of one or more pathways
 * Genes can have zero or more features
 */
class Gene : public GenomicRegion {
public:
	Gene();
	Gene(const char *ensID, uint geneID, const char *chr, uint beg, uint end);
	Gene(const Gene& orig);

	virtual ~Gene();										///< Features are maintained by the chromosome (so we don't delete them here)

	void AddFeature(Feature* feature);				///< Attempts to add a feature to the gene

	std::string EnsemblID();							///< Return the ensemble stable ID

	/**
	 * @brief return all pvalues associated with the gene
	 */
	void CollectPValues(std::vector<float>& pvalues) const;

	/**
	 * @brief Counts signficant features, returning count details
	 */
	int CountSignificantMembers(std::set<Feature*>& features, uint& simpleFeature, uint &complexFeature, uint& sigSimple, uint &sigComplex) const;
	/**
	 * @brief Simply returns the number of significant features found
	 */
	int CountSignificantMembers(std::set<Feature*>& usedFeatures) const;
	/**
	 * @brief Returns number of complex features associated with the gene
	 */
	int CountComplexFeatures();
	/**
	 * @brief Returns the number of simple (1 SNP only) features associated with the gene
	 */
	int CountSimpleFeatures();
	/**
	 * @brief Counts the number of empty features (regions where no datapoints were found from the dataset)
	 */
	int CountEmptyFeatures();

	/**
	 * @brief Generates a detailed report about each feature found within the gene
	 */
	void ComprehensiveReport(std::ostream& os);

	/**
	 * @brief Returns number of features associated with the gene
	 */
	uint FeatureCount() const;

	void CollectSNPs(std::set<uint>& snps);
	/**
	 * Returns the count for each bin associated in the local features
	 */
	void GetFeatureMap(std::multiset<uint>& bins, std::set<uint>& featureIDs);

	/**
	 * @brief Add a group to the gene's group list
	 */
	void AddGroup(uint groupType, uint groupID);	///< Add a group to the gene's group information

	/**
	 * @brief Associate a familiar name with the gene, "i.e. ABCB1"
	 */
	void SetAlias(const char *alias);

	std::set<Feature*> GetFeatures();
	/**
	 * @Brief Return a gene's alias
	 */
	std::string Alias();

	std::multimap<uint, uint> GetGroups();

	void ListFeatures(std::ostream& os);


	void DetailedReport(std::map<uint, uint>& snps, const char *prefix, std::ostream& os);
	void Summarize(std::ostream& os, bool html);
	static bool AllowRedundantFeatures;				///< Allow users to count features that appear in multiple genes (inside same pathway) multiple times
private:
	std::string ensemblID;								///< Ensembl ID
	std::set<Feature*> features;						///< Features associated with this gene
	std::multimap<uint, uint> groups;				///< type -> groupID (type is like KEGG)
	std::string alias;									///< Common name
};





}
#endif	/* _GENE_H */

