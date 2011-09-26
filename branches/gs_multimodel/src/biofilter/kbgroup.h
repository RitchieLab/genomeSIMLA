//
// C++ Implementation: kbgroup.h
//
// Description: Primary functionality associated with a "group" from one of the knowledge sources
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

//
// C++ Interface: group
//
// Description:
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERGROUP_H
#define BIOFILTERGROUP_H

#include "kbentity.h"
#include "kbregion.h"
#include "genegenemodelarchive.h"
#include <soci.h>
#include <soci-sqlite3.h>
namespace Biofilter {
namespace Knowledge {

using namespace soci;

class KbGroup;
typedef std::map<uint, KbGroup*> GroupLookup;

/**
 * @brief Represents a pathway, reaction or other "group" type entity, which can contain regions as well as other groups

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class KbGroup : public KbEntity {
public:
	KbGroup(uint id, uint typeID, bool isDiseaseDependent, const char *name, const char *desc);
    ~KbGroup();


	void Add(KbGroup* child);
	void Add(KbRegion* child);

	/**
	 * @Brief 	Performs a semi-recursive model generation.
	 * @Note  	Once we start looking at 3rd order and greater,
	 * 		  	we will need to reconsider this method
	 * @param other The other "bag" against which we'll be creating new models
	 * @param repo The repository to which we'll write our models to
	 * @return	Returns the number of models associated with this "bag"
	 *			that were added to the repository
	 */
	uint CountChildren() { return children.size();}

	uint GenerateGeneGeneModels(GeneGeneModelArchive& geneArchive, uint pairID, uint maxGeneCount, std::ostream& os);


	/**
	 * @Brief	Produces a list of EVERY unique region (bag).
	 * @note	This is recursively produced
	 */
	void GetAllRegions(std::set<KbRegion*>& allRegions);
	void GetAllGroups(std::set<KbGroup*>& allGroups) const;
	/**
	 * @Brief Ties genes in with a given group
	 * We will build up that list of genes as we go, adding genes as is appropriate
	 * @Return The number of genes that were encountered (including those that were already included)
	 */
	uint AssociateGenes(std::map<uint, KbRegion*>& genes);
	uint AssociateGenes(soci::session& sociDB, std::map<uint,KbRegion*>& genes, SnpManager& snps, std::map<uint, std::string>& aliases, uint popID = 0);
	void AssociateGene(KbRegion* gene);
	/**
	 * @brief Associate snps with the children of this group
	 * @return Returns the number snps that were associated
	 */
	uint AssociateSNPs(SnpManager& snps);

	/**
	 * @Brief ensure that all child groups are loaded and attributed to their parents
	 * @param sociDB the database
	 * @param groups the groups that have already been loaded
	 * @param parents The parents to who this node's work is attributed
	 */
	void LoadAssociations(soci::session& sociDB, std::map<uint, KbGroup*>& groups, std::set<uint> parents);

	/**
	 * @Brief Loads name and description from database
	 */
	void Load(soci::session& sociDB);

	uint ListAssociations(int tabCount, std::ostream& os, uint maxGeneCount);
	uint GraphAssociations(std::ostream& os, uint maxGeneCount);

	//void Release();

	void MarkProcessed(bool isProcessed);
	void AddGroupTypes(std::set<uint>& types);


	typedef std::set<KbGroup*>::iterator iterator;
	iterator begin();
	iterator end();
	
	std::set<KbRegion*> GetRegions();


	enum DiseaseDependentRelationshipType {
		AllModels,
		DD_GroupLevel,
		DD_Only
	};
	static DiseaseDependentRelationshipType DiseaseDependentRelationship;
	static bool CollapseAssociationReport;
protected:
	/**
 	 * If there is the possibility of instantiating more than one of each type, we probably will need
	 * to use our own comparison operator to dereference the pointer...but, I would hope we never have
	 * more than one instance of each item.
	 */
	std::set<KbGroup*> children;			///< Child groups
	std::set<KbRegion*> regions;			///< Regions associated with local node
	uint groupType;
	bool isDiseaseDependent;				///< Distinguish between DI/DD groups
	//std::set<uint> groupType;				///< This is used for reporting the implication index
	bool genesAssociated;			///< Short circuit the loading process for groups that appear more than once in a tree
};

inline
std::set<KbGroup*>::iterator KbGroup::begin() {
	return children.begin();
}

inline
std::set<KbGroup*>::iterator KbGroup::end() {
	return children.end();
}
	
inline
std::set<KbRegion*> KbGroup::GetRegions() {
	return regions;
}

}
}
#endif
