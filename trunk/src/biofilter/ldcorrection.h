//
// C++ Interface: ldcorrection
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BIOFILTERLDCORRECTION_H
#define BIOFILTERLDCORRECTION_H
#include <map>
#include <string>
#include "kbregion.h"
#include "regionspline.h"
#include "utility/rbtree.h"

namespace Biofilter {

typedef Utility::RBTree<uint, RegionSpline*> RegionBoundaries;
typedef Utility::RBTreeNode<uint, RegionSpline*> RegionBoundariesNode;
typedef std::map<uint, SNP_Details*> SnpLookup;

/**
@Brief Adds references to extended region boundaries based on specified LD values

	@author Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>
*/
class LdCorrection : public SnpManager {
public:
    LdCorrection();
    ~LdCorrection();

	/**
	 * @Brief Update the region bounds 
	 */
	void Process(soci::session& sociDB, const char *variationFilename);

	/**
	 * @brief Load the various LD settings from file, cfg
	 */
	void LoadConfiguration(soci::session& sociDB, const char *cfg);

	/**
	 * @Brief Initialize snps based on their presence in the set, snps
	 * @PARAM chrom - Which chromosome
	 * @PARAM snps  - The set containing rs IDs to be loaded
	 * @PARAM variation_filename - source file containing snp positions
	 */
	uint InitSNPs(int chrom, std::set<uint>& snps, const char *variation_filename);

	//void PrintSNPs(ostream& os);

protected:
	void LoadGenes(soci::session& sociDB, int chrom);
	void PurgeRegions();
	/**
	 * @Brief Initialize the population, based on the default values
	 */
	int InitPopulation(soci::session& sociDB, const char *pop, const char *popDesc);

	/**
	 * @brief Process the LD from file, filename
	 */
	void ProcessLD(soci::session& sociDB, int chrom, const char *ldSource, const char *variationFilename);

	std::string GetPopDescription(const char *stat, float threshold);
	std::string GetPopulationName(const char *type, float threshold);
	int GetPopID(soci::session& sociDB, const char *type, float threshold);
	/**
	 * @Brief Load ld boundaries from a stream (into the vector)
	 */
	void LoadValuesDP(soci::session& sociDB, std::istream& os);
	void LoadValuesRS(soci::session& sociDB, std::istream& os);


	std::string ldSourceFile;			///< Filename used to load LD

	// The following will be added to the database 
	std::string ldName;					///< Short name used to symbolically refer to ld settings 
	std::string ldComment;				///< Comment describing details about this population
	int popID;							///< The population ID associated with the current population

	std::vector<std::string> ldFilenames;			///< The files to be used to load LD
	RegionBoundaries regions;			///< Region Splines associated with the current ld file

	//SnpLookup snpData;					///< rs->SNP
	std::string snpDataFilename;				///< The name of the variation data

	//SNP_Details *first;							///< first SNP in the snpData structure (according position)
	//SNP_Details *last;							///< the last SNP in the snpData
};




}

#endif
