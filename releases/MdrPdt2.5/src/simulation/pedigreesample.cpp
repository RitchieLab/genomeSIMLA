//
// C++ Implementation: pedigreesample
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreesample.h"
#include "poolmanager.h"
#include <iomanip>

namespace Simulation {

using namespace Utility;

PedigreeSample::PedigreeSample(uint aff, uint unaff, float genoError, float phenoError, float missing, uint extras) : 
		Sample(genoError, phenoError, missing), affectedSibCount(aff), unaffectedSibCount(unaff), extraSibs(extras), pedID(0) { }


PedigreeSample::~PedigreeSample() { }

void PedigreeSample::ResetPedID() { 
	pedID=1;
}

void PedigreeSample::DrawIndividual(PoolManager &pools, Individual &person) {
	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();
	while (pool) {
		pool->DrawIndividual( &person);
		pool  = itr.GetNext();
	}
}

void PedigreeSample::ReturnIndividual(PoolManager &pools, Individual &person) {
	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();
	while (pool) {
		pool->ReturnIndividual( &person);
		pool  = itr.GetNext();
	}

}	

void PedigreeSample::BuildSample(PoolManager &pools, ModelManager& models, uint familyCount) {
	people.clear();

	PenetranceModel *model;
	uint familiesRemaining = familyCount;
	uint modelCount = models.GetModelCount();
	assert(modelCount > 0);

	for (uint i=0; i<modelCount - 1; i++) {
		model = models.GetModel(i);
		uint percentage = (uint)((float)familyCount * model->GetProbability());
		familiesRemaining-=percentage;
		BuildSample(pools, model, percentage);
	}

	//OK, let's just finish of the last bunch with the 
	model = models.GetModel(modelCount - 1);
	BuildSample(pools, model, familiesRemaining);
	
}

void PedigreeSample::BuildSample(PoolManager &pools, PenetranceModel* model, uint individualCount) {
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models
	uint childCount = affectedSibCount + unaffectedSibCount + extraSibs;
	uint chromCount = pools.Size();
	
	Purge();
	
	familyStartPositions.clear();
	while (individualCount-- > 0) {
		pedID++;
		bool validPedigree=false;

		uint attempts = 0;

		while (!validPedigree) {
			if (attempts++ > 10000) {
				cout<<"Unable to find a complete dataset after 10000 tries\n";
				attempts=0;
			}
			uint id=1;

			uint affecteds = 0;

			Individual *dad = new Individual(id++, pedID, chromCount);
			Individual *mom = new Individual(id++, pedID, chromCount);

			Individual *children[childCount];
			//for (uint i=0; i<childCount; i++) {
			//	children[i]=new Individual(id, pedID, chromCount);
			//	children[i]->SetPedigreeMeta(2, 1);
			//}

			DrawIndividual(pools, *mom);
			DrawIndividual(pools, *dad);
			//Build up the children
			for (uint i=0; i<childCount; i++) {
				//children[i]->Init(id++, pedID, chromCount);
				//This is broken! I am not using crossover to generate children!
				children[i]=dad->Cross(mom, id++);
				
				if (children[i]->ApplyStatus(model)) 
					affecteds++;
			}

			//Test for appropriate familial shape
			if (affecteds >= affectedSibCount && childCount - affecteds >= unaffectedSibCount) {
				validPedigree=true;
				//Record the start position in our special array
				familyStartPositions.push_back(people.size());
				people.push_back(dad);
				people.push_back(mom);
				uint unaffecteds=unaffectedSibCount;
				affecteds=affectedSibCount;
				for (uint i=0; i<childCount; i++) {	
					if (children[i]->IsAffected()) {
						if (affecteds>0) {
							people.push_back(children[i]);
							affecteds--;
						}
					}
					else {
						if (unaffecteds>0) {
							people.push_back(children[i]);
							unaffecteds--;
						}
					}
				}
			}
			else {
				ReturnIndividual(pools, *dad);
				delete dad;
				ReturnIndividual(pools, *mom);
				delete mom;
				for (uint i=0; i<childCount; i++) 
					delete children[i];

			}
		}
	}

}

/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 */
int PedigreeSample::WriteDataset(ostream &os, uint *gtCount) {
	vector<uint> startPoints = familyStartPositions;
	random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	
	uint count = startPoints.size();
	uint familySize = 2 + affectedSibCount + unaffectedSibCount;



	if (familySize == 2)
		familySize +=extraSibs;
	uint pos = 0;
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		pos = startPoints[i];
		for (uint o=0;o<familySize; o++) {
			Individual *person=people[pos++];
			if (person->DoIncludeInDataset()) {
				observations++;
				person->WritePedigree(os, gtCount);
			}
		}
	}
	return observations;

}

void PedigreeSample::WritePhased(ostream &os) {
	uint count = familyStartPositions.size();

	uint idx = 1;
	for (uint i=0; i<count; i++) {
		people[familyStartPositions[i]]->WritePhased(os, idx++);
		people[familyStartPositions[i] + 1]->WritePhased(os, idx++);
	}
}

/** I need to figure out how to id the mom and dad */
void PedigreeSample::ApplyPhenocopyError(PoolManager *pools, uint familyCount) {
	uint familySize = 2 + affectedSibCount + unaffectedSibCount;

	uint sampleSize = people.size();
	uint errorCount = (uint)((float)(familyCount * affectedSibCount) * phenoError);

	BitSetType tested(sampleSize, false);

	Individual *individual=NULL;
	for (uint i=0; i<errorCount; i++) {
		bool isAffected = false;
		uint idx;
		while (!isAffected) {
			idx = Utility::Random::globalGenerator((int)sampleSize);
			if (!tested[idx]) {
				individual = people[idx];
				isAffected = individual->IsAffected();
				tested[idx]=true;
			}
		}

		uint dad = idx - (idx % familySize);
		uint mom = dad+1;
		uint newID = individual->GetID() + 100;
		delete individual;
		Individual *d = people[dad];
		Individual *m = people[mom];

		individual = d->Cross(m, newID);
		individual->SetStatus( true );

		people[idx]=individual;
	}
}



void PedigreeSample::GenerateReport(ostream &os, uint padding) {
	os<<setw(padding - 5)<<"Pedigree Sample : "<<affectedSibCount<<"A/"<<unaffectedSibCount<<"U "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(genoError*100.0)<<"% Genotype Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(phenoError*100.0)<<"% Phenocopy Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(missingData*100.0)<<"% Missing Data "<<endl;
}
















void PedigreeMixedSample::BuildSample(PoolManager &pools, PenetranceModel* model, uint individualCount) {
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models
	uint chromCount = pools.Size();
	
	Purge();
	
	familyStartPositions.clear();
	while (individualCount-- > 0) {
		pedID++;
		bool validPedigree=false;

		uint attempts = 0;

		while (!validPedigree) {

			//Random number generator returns a value from 0 to n-1
			uint childCount = Utility::Random::globalGenerator((int)extraSibs) + 1;

			if (attempts++ > 10000) {
				cout<<"Unable to find a complete dataset after 10000 tries\n";
				attempts=0;
			}
			uint id=1;

			uint affecteds = 0;

			Individual *dad = new Individual(id++, pedID, chromCount);
			Individual *mom = new Individual(id++, pedID, chromCount);

			Individual *children[extraSibs];

			DrawIndividual(pools, *mom);
			DrawIndividual(pools, *dad);
			//Build up the children
			for (uint i=0; i<extraSibs; i++) {
				//children[i]->Init(id++, pedID, chromCount);
				//This is broken! I am not using crossover to generate children!
				children[i]=dad->Cross(mom, id++);
				
				if (i < childCount) {
					if (children[i]->ApplyStatus(model))  {
						affecteds++;
					}
				}
				else 		//This will let us skip this child
					children[i]->DoIncludeInDataset(false);
			}

			//Test for appropriate familial shape
			if (affecteds > 0) {
				validPedigree=true;
				//Record the start position in our special array
				familyStartPositions.push_back(people.size());
				people.push_back(dad);
				people.push_back(mom);
				for (uint i=0; i<extraSibs; i++) {	
					people.push_back(children[i]);
				}
			}
			else {
				ReturnIndividual(pools, *dad);
				delete dad;
				ReturnIndividual(pools, *mom);
				delete mom;
				for (uint i=0; i<extraSibs; i++) 
					delete children[i];

			}
		}
	}

}


}
