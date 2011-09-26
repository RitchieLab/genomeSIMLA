//
// C++ Interface: kbmetagroup.h
//
// Description: Primary functionality associated with a single biofilter knowledge source
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef BIOFILTERMETAGROUP_H
#define BIOFILTERMETAGROUP_H
#include "kbentity.h"
#include "genegenemodelarchive.h"
#include "kbgroup.h"
namespace Biofilter {
namespace Knowledge {
/**
	@brief Organizes groups by their primary group type (such as Reactome, GO, etc)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class KbMetaGroup : public KbEntity{
public:
	typedef std::map<uint, KbGroup*>::iterator iterator;
    KbMetaGroup(uint id,  bool isDiseaseDependent, const char *name, const char *date);
	KbMetaGroup(uint id, bool isDiseaseDependent = false);
    ~KbMetaGroup();


	/**
	 * @Brief Loads the groups from the database
	 * @Param sociDB the database session
	 * @Param maxComponentThreshold Groups larger than this will be removed from memory
	 * @param inclusions if not empty, only IDs in this list will be used (this is a comma separated list)
	 */
	uint LoadGroups(soci::session& sociDB, int maxComponentThreshold, const char *inclusions, GroupLookup& grpLookup);
	uint LoadGroups(soci::session& sociDB, int maxThesh, std::vector<uint>& roots, GroupLookup& grpLookup);
	/**
	 * @Brief Ties genes in with a given group
	 * We will build up that list of genes as we go, adding genes as is appropriate
	 * @Return The number of genes that were encountered (including those that were already included)
	 */
	uint AssociateGenes(soci::session& sociDB, std::map<uint,KbRegion*>& genes, SnpManager& snps, std::map<uint, std::string>& aliases, uint popID = 0);
	void AddGroup(KbGroup *group);

	void ListAssociations(std::ostream& os, uint maxGeneCount);
	void GraphAssociations(std::ostream& os, uint maxGeneCount);

	uint GenerateGeneGeneModels(GeneGeneModelArchive& archive, int maxGeneCount, std::ostream& os);

	iterator begin();
	iterator end();

	uint GetGroupCount();

	bool IsDiseaseDependent();
protected:
	std::string date;				///< This will be available if we are reporting
	GroupLookup groups;	///< This needs to be a lookup type
	bool isDiseaseDependent;				///< Distinguish between DI/DD groups

};



}
}

#endif
