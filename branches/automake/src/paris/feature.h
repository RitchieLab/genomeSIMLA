/* 
 * File:   feature.h
 * Author: torstees
 *
 * Features represent genomic regions within genes.
 *		Comprised of one or more points from a dataset (SNPs)
 *
 * Created on January 5, 2010, 11:52 AM
 */

#ifndef PARIS_FEATURE_H
#define	PARIS_FEATURE_H

#include <set>
#include <map>
#include <string>
#include "genomicregion.h"
#include <soci.h>


namespace soci {
namespace details {

template <>
struct exchange_traits<unsigned int>
{
        typedef basic_type_tag type_family;
        enum { x_type = x_integer };
};
}
}


namespace Paris {



/**
 * Features aggregate the signals which are determined to be synonymous due to hi LD correlation
 */
class Feature : public GenomicRegion {
public:
	Feature();
	Feature(uint id, const char *chr, uint beg, uint end);
	Feature(const Feature& orig);
	virtual ~Feature();
	
	/**
	 *	@brief Attempts to add a value to the feature
	 * @return true if the snp does fit within the boundaries
	 */
	bool AddValue(uint snpIndex, const char *chrom, uint position, float pvalue);

	/**
	 * @brief returns the number of significant members found within the local feature
	 */
	int CountSignificantMembers() const;	///< Returns sigCount (calculating it, if it hasn't been done yet)

	/**
	 * @Brief Writes details of the feature to the stream
	 */
	void DetailedReport(std::map<uint, uint>& snps, const char *prefix, std::ostream& os, uint& totalSig, uint &totalNSig);	///< Report all SNPs associated with the feature

	void WriteBinReport(std::ostream& os);

	/**
	 * @brief SNP IDX (RS) -> score
	 */
	std::map<uint, float> GetPValues();

	std::set<uint> GeneIDs();					///< Return the set of gene IDs
	uint BinIndex();								///< Return the bin Index
	void BinIndex(uint index);					///< Allow the bin to set the index
	uint FeatureSize() const;					///< Return number of SNPs associated with the feature
	static bool IgnorePValueOfZero;			///< Allow the user to ignore pvalues of zero (or not). These will be counted as not significant
private:
	std::set<uint> geneIDs;						///< List of gene IDs which this feature is associated
	uint binIndex;									///< Which bin will this feature be found
	std::map<uint, float> pscores;			///< SNP Idx -> score
	uint sigCount;									///< Count of significant values
};

/**
 * This can be used for sorting in the initial set before binning begins. Once
 * we have a sorted vector or set based on this functor, we can just iterate from
 * beginning to end in order to build up each bin
 */
struct SortByFeatureSize {
	bool operator()(const Feature* left, const Feature* right) {
		if (left->FeatureSize() == right->FeatureSize())
			return left < right;
		return left->FeatureSize() < right->FeatureSize();
	}
};




}


#endif	/* _FEATURE_H */

