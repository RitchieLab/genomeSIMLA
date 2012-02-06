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
namespace StatusModel {
using namespace SimPen;


bool 	GAModel::doRunGA = true;						///<Do we want to use GA to find the best model?
int  	GAModel::seed = 1347;							///<the seed to be used for the ga
int  	GAModel::tries = 5;    							///<How many attempts do we want to do before we give up?
float 	GAModel::fitnessThreshold = 0.95;				///<Threshold before retrying the GA search



GAModel::GAModel()  { }


GAModel::~GAModel(){

}

string GAModel::GetModelConfiguration() {
	stringstream ss;
	ss<<"DEFINE_MODEL SIMPEN LABEL \""<<filename<<"\" ";
	ModelLociArray::iterator itr = loci.begin();
	ModelLociArray::iterator end = loci.end();

	while (itr != end) {
		ss<<itr->label<<" ";
		itr++;
	}
	return ss.str();
}

void GAModel::GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci) {
	os<<"<P><H3>SIMPEN Based Model (TBD):</H3>\n";
	if (modelSize != diseaseLoci.size()) 
		os<<"<B>Invalid Locus List Specified</B>";
	ReportDiseaseLoci(os, diseaseLoci);

	ReportPenetranceTable(os);

}

void GAModel::Load() {
	//let's make a local copy so we can iterate over seeds when we figure out what is successful or not
	int seed = this->seed;
	float score =0.0;

	float bestFitness = 0.0;

	string modelFilename 	= filename + ".1";
	string bestModel 		= modelFilename;
	

	if (doRunGA) {		
		score = 0.0;
		
		

		/**
		 * @todo I need to break up the model selection so that we can actually save/write the best
		 * score over the iterations. Right now, it will write out the last one that was produced.
	 	 */
		cout<<"Searching for a suitable penetrance table. \n";
		cout<<"* The best model will be selected from "<<tries<<" or the first to exceed "<<fitnessThreshold<<"\n";

		for (int n=0; n<tries && score<fitnessThreshold; n++) {
			score = run_simpen(filename.c_str(), modelFilename.c_str(), loci, seed++, true);	
			if (score> bestFitness) {
				bestFitness=score;
				bestModel=modelFilename;
			}				
			else if (score == 0.0) {
				cout<<"!! There was a serious problem when trying to load your model configuration file: "<<filename<<"\n";
				cout<<"Please verify that the file exists and is has no errors and try again.\n";
				return;
			}
			cout<<"\tFitness - "<<score<<" ("<<modelFilename + string(".1.smod")<<") -- Seed: "<<seed-1<<"\n";
			char b[12];
			sprintf(b, "%d", n+1);
			modelFilename = filename + "." + b;
		}
	} 
	if (bestFitness > 0.0) {
		if (bestFitness < 1.0) 
			cout<<"*************** The best model had a fitness below 1.0 (it didn't meet all of the criterion asked for)\n";
		
		string penFile = bestModel + string(".1.smod");
		PenFileModel::Load(penFile.c_str());
	}
}

}

}
