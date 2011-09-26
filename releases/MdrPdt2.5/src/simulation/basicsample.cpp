//
// C++ Implementation: population
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "utility/types.h"
#include "basicsample.h"
#include "poolmanager.h"
#include <iomanip>

namespace Simulation {

using namespace std;
using namespace Utility;

void BasicSample::BuildSample(PoolManager &pools, PenetranceModel* model, uint individualCount) {
	uint affecteds = (uint)((float)individualCount * percentAffected);
	uint unaffecteds = individualCount - affecteds;
	
	uint id=0;					///<We have to keep up with individuals separately
	//For now, lets make this as simple as possible. However, eventually, we will need to consider balancing
	//the different models
	while (affecteds > 0 || unaffecteds > 0) {
		id++;
		bool validIndividual=false;

		while (!validIndividual) {
			Individual *newPerson = new Individual(id, id);
			PoolManager::Iterator itr = pools.GetIterator();
			ChromPool *pool = itr.GetNext();
			while (pool) {
				pool->DrawIndividual( newPerson);
				pool=itr.GetNext();
			}
			
			//At this point, we have all the data associated with our person
			//Time to assign status 
			if (newPerson->ApplyStatus(model)) {
				if (affecteds > 0) {
					affecteds--;
					people.push_back(newPerson);
					validIndividual = true;
				}
			}
			else {
				if (unaffecteds > 0) {
					unaffecteds--;
					people.push_back(newPerson);
					validIndividual = true;
				}
			}
		}
	}


}
void Sample::ApplyMissingData(PoolManager *mgr, uint familycount) {
	uint sampleSize = people.size();

	if (sampleSize < 0)
		return;
	PoolManager::Iterator itr = mgr->GetIterator();

	ChromPool *ch = itr.GetNext();
	//Iterate over each of the chromosomes
	while (ch) {
		//uint chrCount = ch->GetID();

		uint locusSize = ch->GetLociCount();

		//Iterate through each locus and modify genoError% of individuals
		for (uint locus=0; locus<locusSize; locus++) { 

			uint errorCount = (uint)((float)sampleSize * missingData);
			Utility::BitSetType alreadyUsed(sampleSize, false);

			while (errorCount > 0) {
				int idx = Utility::Random::globalGenerator((int)sampleSize);
				if (!alreadyUsed[idx]) {
					people[idx]->ClearLocus(ch->GetID(), locus);
					alreadyUsed[idx]=true;
					errorCount--;
				}
			}
		}
		ch = itr.GetNext();
	}
}


void Sample::ApplyErrors(PoolManager *mgr, uint familycount) {
	ApplyGenocopyError( mgr, familycount );
	ApplyPhenocopyError( mgr, familycount );
	ApplyMissingData( mgr, familycount );
}

void BasicSample::ApplyPhenocopyError(PoolManager *pools, uint familycount) {
	uint sampleSize = people.size();
	uint errorCount = (uint)((float)sampleSize * phenoError);
	uint topID = people[sampleSize - 1]->GetID();

	BitSetType tested(sampleSize, false);

	Individual *individual=NULL;
	for (uint i=0; i<errorCount; i++) {
		bool isAffected = false;
		int idx;
		while (!isAffected) {
			idx = Utility::Random::globalGenerator((int)sampleSize);
			if (!tested[idx]) {
				individual = people[idx];
				isAffected = individual->IsAffected();
				tested[idx]=true;
			}
		}
				
		delete individual;
		individual = new Individual(++topID, topID);

		PoolManager::Iterator itr = pools->GetIterator();
		ChromPool *pool = itr.GetNext();
		while (pool) {
			pool->DrawIndividual( individual);
			pool=itr.GetNext();
		}
		individual->SetStatus( true );
	}
}


/**
 * Apply the error according to the locus rate
 */
void Sample::ApplyGenocopyError(PoolManager *mgr, uint familycount) {
	uint sampleSize = people.size();

	if (sampleSize < 0)
		return;
	PoolManager::Iterator itr = mgr->GetIterator();

	ChromPool *ch = itr.GetNext();
	//Iterate over each of the chromosomes
	while (ch) {
//		uint chrCount = ch->GetID();

//		for (uint chrID=0; chrID<chrCount; chrID++) {
			uint locusSize = ch->GetLociCount();

			//Iterate through each locus and modify genoError% of individuals
			for (uint locus=0; locus<locusSize; locus++) { 
				//float genoError = ch->GetErrorRate(locus);
				uint errorCount = (uint)((float)sampleSize * genoError);
				Utility::BitSetType alreadyUsed(sampleSize, false);
	
				while (errorCount > 0) {
					uint idx = Utility::Random::globalGenerator((int)sampleSize);
					if (!alreadyUsed[idx]) {
						people[idx]->ChangeGenotype(ch->GetID(), locus, 0);
						alreadyUsed[idx]=true;
						errorCount--;
					}
				}
			}
//		}
		ch = itr.GetNext();
	}
}

/**
void Sample::ApplyGenocopyError(float genoError, PoolManager *mgr) {
	assert(genoError < 1.0);
	uint sampleSize = people.size();

	if (sampleSize < 0)
		return;
	PoolManager::Iterator itr = mgr->GetIterator();

	ChromPool *ch = itr.GetNext();
	//Iterate over each of the chromosomes
	while (ch) {
		uint chrCount = ch->GetID();

		for (uint chrID=0; chrID<chrCount; chrID++) {
			uint locusSize = ch->GetLociCount();

			//Iterate through each locus and modify genoError% of individuals
			for (uint locus=0; locus<locusSize; locus++) { 
				uint errorCount = (uint)((float)sampleSize * genoError);
				int errDir = ch->GetErrorDir(locus);
			
				Utility::BitSetType alreadyUsed(sampleSize, false);
	
				while (errorCount > 0) {
					uint idx = Utility::Random::globalGenerator((long)sampleSize);
					if (!alreadyUsed[idx] && people[idx]->ChangeGenotype(chrID, locus, errDir))
						alreadyUsed[idx]=true;
				}
			}
		}
		ch = itr.GetNext();

	}		
} */

void Sample::Purge() {
	for (uint i=0; i<people.size(); i++)
		delete people[i];
	people.clear();
}

void BasicSample::BuildSample(PoolManager &pools, ModelManager& models, uint individualCount) {
	Purge();

	PenetranceModel *model;
	uint individualsRemaining = individualCount;
	uint modelCount = models.GetModelCount();
	assert(modelCount > 0);

	for (uint i=0; i<modelCount - 1; i++) {
		model = models.GetModel(i);
		uint percentage = (uint)((float)individualCount * model->GetProbability());
		individualsRemaining-=percentage;
		BuildSample(pools, model, percentage);
	}

	//OK, let's just finish of the last bunch with the 
	model = models.GetModel(modelCount - 1);
	BuildSample(pools, model, individualsRemaining);
	
}


/**
 * @brief Dump the contents of the sample population to the stream
 * @param os The stream to which the sample will be written
 */
int BasicSample::WriteDataset(ostream &os, uint *gtCounts) {
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WriteMDR(os, gtCounts);

	return count;
}

/**
 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
 * @param os The stream to which the sample will be written	
 */
void BasicSample::WritePhased(ostream &os) {
	uint count = people.size();
	for (uint i=0; i<count; i++) 
		people[i]->WritePhased(os, i+1);

}
void BasicSample::GenerateReport(ostream &os, uint padding) {
	os<<setw(padding - 5)<<"Case Control Sample : "<<setprecision(2)<<(100.0 * percentAffected)<<"% "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(genoError*100.0)<<"% Genotype Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(phenoError*100.0)<<"% Phenocopy Error "<<endl;
	os<<setw(padding - 15)<<""<<setprecision(2)<<setw(8)<<(missingData*100.0)<<"% Missing Data "<<endl;
}



}
