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
#include "penfilemodel.h"
#include "poolmanager.h"
#include "utility/types.h"

namespace Simulation {

namespace StatusModel {
using namespace std;

/**
@brief Adaptor for generating a penetrance table using the simpen library. 
Has ability to use the GA to produce the penetrance table for status based on a specialized configuration file and the allele frequencies from the actual data

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GAModel : public PenFileModel {
public:
	GAModel();
	~GAModel();

	/**
	 * @brief Load the model at the local file, filename
	 */
	void Load();
	string GetModelConfiguration();

//	virtual bool Init(istream &i, PoolManager* pools);

	string Details();
	void GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci);
	static bool doRunGA;					///<Do we run the GA or skip to the file used
	static int seed;						///<the seed to be used for the ga
	static int tries;						///<How many attempts do we want to do before we give up?
	static float fitnessThreshold;			///<How good of a fitness do we require before starting again
	int GetType() { return 4; }
protected:
	vector<Utility::LocusID> locusIDs;

};

inline
string GAModel::Details() {
	return string("Simpen Model ") + filename;
}

}

}

#endif
