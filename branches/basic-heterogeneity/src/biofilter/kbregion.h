//
// C++ Interface: region
//
// Description:
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERREGION_H
#define BIOFILTERREGION_H

#include "kbentity.h"
#include "region.h"
#include "utility/strings.h"
#include "snpmanager.h"
namespace Biofilter {
namespace Knowledge {
/**
	@Brief Represents a single genomic region (such as a gene or linkage region)
	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class KbRegion : public KbEntity, public Region {
public:
	struct SnpType {
		typedef enum Type {
			Interior,			///< Within true gene boundaries
			Flanking,			///< Outside actual boundaries, but within extended boundaries
			Exterior,			///< Outide even extended boundaries
			Unknown				///< A problem occured and we don't know. Maybe this is the result of trying to report on more than one SNP at a time?
		} Type;
	};

	KbRegion();
    KbRegion(uint geneID, uint start, uint end, const char *chrom, const char *ensembl, const char *desc, SnpManager* snpMgr);
	KbRegion(const KbRegion& other);

    ~KbRegion();

	/**
	 * @brief consider adding a SNP to the snp bag
	 * @return T/f indicating success of adding the snp
	 */
	bool AddSnp(uint snp);

	void SetSnpManager(SnpManager* snpMgr);

	void GetBounds(uint& start, uint& end);
	void GetOrigBounds(uint& start, uint& end);

	void SetBounds(uint start, uint end);
	void SetOrigBounds(uint start, uint end);
	/**
	 * @Brief Associate any SNPs with the local region based on location
	 */
	uint AssociateSNPs();

	uint ListAssociations(int tabCount, std::ostream& os);

	uint SnpCount();

	int GetSnpCoverage(std::set<uint>& snpList, std::set<SNP_Details>& snpDetails);
	uint CollectSnpDetails(std::set<SNP_Details>& collection);

	/**
	 * @brief Counts the number of RS Ids found int he set, snpList
	 */
	uint SnpCount(std::set<uint>& snpList);

	SnpType::Type GetSnpTypeByRS(uint snpPosition);
	/**
	 * @Brief Returns the matching SNPs found in the snps and the local snp list
	 */
	std::set<uint> GetSnpCoverage(std::set<uint>& snps);

	uint Start();
	uint OriginalStart();
	uint End();
	uint OriginalEnd();
	std::string Chromosome();
	std::string RegionName();

	std::set<uint> SNPs();
protected:
	uint startPosition;				///< Position representing the start of the region
	uint endPosition;					///< Position representing the end
	uint origStart;					///< This is the official start position (not based on LD extension)
	uint origEnd;						///< This si the official end position
	std::string chromosome;			///< Chromosome (used in reporting)
	SnpManager *snpLookup;			///< Used for translating back to rs IDs
};



}

}

#endif
