//
// C++ Interface: gamodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONGAMODEL_H
#define SIMULATIONGAMODEL_H

#include <string>
#include "penetrancemodel.h"
#include "poolmanager.h"

namespace Simulation {

using namespace std;

/**
@brief Adaptor for generating a penetrance table using the simpen library. 
Has ability to use the GA to produce the penetrance table for status based on a specialized configuration file and the allele frequencies from the actual data

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GAModel : public PenetranceModel {
public:
	GAModel(uint modelID, const char *name, float prob, vector<int> chr, vector<int> loci);
	~GAModel();

	/**
	 * @brief Load the model at the local file, filename
	 */
	void Load();
	
	static bool doRunGA;					///<Do we run the GA or skip to the file used
	static int seed;						///<the seed to be used for the ga
	static int tries;						///<How many attempts do we want to do before we give up?
	static float fitnessThreshold;			///<How good of a fitness do we require before starting again

	void SetPoolManager(PoolManager *mgr);	///<This is required- we need to know who to talk to for chromosome data
protected:
	string gaConfigFilename;				///<Where we expect to found the configuration details
	//string modelFilename;					///<Where we will store the penetrance table generated
	int *loci;								///<The disease loci
	int *chr;								///<The chromosome where each loci resides
	int locCount;							///<The number of loci in the final model
	PoolManager *poolMgr;					///<This is just a copy- so no need to delete it 
};

inline
void GAModel::SetPoolManager(PoolManager *mgr) {
	poolMgr = mgr;
}

}


#endif
