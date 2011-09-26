//
// C++ Implementation: gamodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gamodel.h"
#include <simpen/simpen.h>

namespace Simulation {

using namespace SimPen;


bool 	GAModel::doRunGA = true;						///<Do we want to use GA to find the best model?
int  	GAModel::seed = 1347;							///<the seed to be used for the ga
int  	GAModel::tries = 5;    							///<How many attempts do we want to do before we give up?
float 	GAModel::fitnessThreshold = 0.95;				///<Threshold before retrying the GA search


GAModel::GAModel(uint modelID, const char *name, float prob, vector<int> chr, vector<int> loci) : 
		PenetranceModel(modelID, prob), locCount(loci.size()), poolMgr(NULL) {
	//modelFilename = string(name) + ".model";
	gaConfigFilename = string(name);
		
	this->loci = new int[locCount];
	this->chr  = new int[locCount];

	cout<<"Building Disease Loci for Model "<<modelID<<":\n\t";
	for (int i=0; i<locCount; i++) {
		if (i>0)
			cout<<"x";
		cout<<chr[i]<<":"<<loci[i];
		this->loci[i]=loci[i];
		this->chr[i] =chr[i];
	}
	cout<<"\n";
}


GAModel::~GAModel(){
	if (loci)
		delete[] loci;
	if (chr)
		delete[] chr;
}


void GAModel::Load() {
	//Assertions here mean you didn't first set the pool manager before loading
	assert(poolMgr != NULL);
	
	//let's make a local copy so we can iterate over seeds when we figure out what is successful or not
	int seed = this->seed;
	float score =0.0;

	float bestFitness = 0.0;

	string modelFilename 	= gaConfigFilename + ".1";
	string bestModel 		= modelFilename;
	

	if (doRunGA) {		
		score = 0.0;
		SimPen::LocusArray l;
		
		//Set up the model frequencies
		for (int i=0; i<locCount; i++){ 
			float af1, af2;
			poolMgr->GetAlleleFrequency(chr[i], loci[i], af1, af2);
			l.push_back(DiseaseLocus(chr[i], loci[i], af1, af2));
			cout<<"\tSNP : "<<chr[i]+1<<":"<<loci[i]+1<<" "<<af1<<", "<<af2<<"\n";
		}
		

		/**
		 * @todo I need to break up the model selection so that we can actually save/write the best
		 * score over the iterations. Right now, it will write out the last one that was produced.
	 	 */
		cout<<"Searching for a suitable penetrance table. \n";
		cout<<"* The best model will be selected from "<<tries<<" or the first to exceed "<<fitnessThreshold<<"\n";

		for (int n=0; n<tries && score<fitnessThreshold; n++) {
			score = run_simpen(gaConfigFilename.c_str(), modelFilename.c_str(), l, seed++, true);	
			if (score> bestFitness) {
				bestFitness=score;
				bestModel=modelFilename;
			}				
			else if (score == 0.0) {
				cout<<"!! There was a serious problem when trying to load your model configuration file: "<<gaConfigFilename<<"\n";
				cout<<"Please verify that the file exists and is has no errors and try again.\n";
				return;
			}
			cout<<"\tFitness - "<<score<<" ("<<modelFilename + string(".1.smod")<<") -- Seed: "<<seed-1<<"\n";
			char b[12];
			sprintf(b, "%d", n+1);
			modelFilename = gaConfigFilename + "." + b;
		}
	} 
	if (bestFitness > 0.0) {
		if (bestFitness < 1.0) 
			cout<<"*************** The best model had a fitness below 1.0 (it didn't meet all of the criterion asked for)\n";
		
		string penFile = bestModel + string(".1.smod");
		PenetranceModel::Load(penFile.c_str());
	}
}

}
