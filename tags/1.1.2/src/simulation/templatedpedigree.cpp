//
// C++ Implementation: templatedpedigree
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "templatedpedigree.h"
using namespace std;

namespace Simulation {
namespace PedigreeTemplates {

float 	TemplatedPedigree::tolerance 	= 1.0;
int 	TemplatedPedigree::maxAttempts 	= 1000;

TemplatedPedigree::TemplatedPedigree() : totalAttempts(0) { }

TemplatedPedigree::~TemplatedPedigree() { 
	PurgeIndividuals(founders);
	IndividualLookup::iterator itr = individuals.begin();
	IndividualLookup::iterator end = individuals.end();
	while (itr != end) {
		delete itr++->second;
	}
	individuals.clear();
	PurgeMembers();
}

int TemplatedPedigree::Load(const char *filename) {
	ifstream file(filename);
	if (!file.good()) {
		cerr<<"Unable to open pedigree reference file: "<<filename<<". Unable to continue.";
		exit(1);
	}

	while (file.good()) {
		char line[4096];
		file.getline(line, 4096);
		TemplatedIndividual *ind = new TemplatedIndividual();
		if (ind->Load(line)) {
			if (individuals.find(ind->GetID()) != individuals.end())
				delete individuals[ind->GetID()];
			individuals[ind->GetID()] = ind;
			origOrder.push(ind->GetID());
			statusCounts.Append(ind->status);
		}
		else
			delete ind;
	}
	cout<<"Template Pedigree Loaded: "<<filename<<".\n";
	statusCounts.Header(cout);
	statusCounts.Report(cout, "Template Layout");
	return individuals.size();
}

void TemplatedPedigree::PurgeIndividuals(std::vector<Individual*>& individuals) {
	vector<Individual*>::iterator itr = individuals.begin();
	vector<Individual*>::iterator end = individuals.end();
	while (itr != end) 
		delete *itr++;
	individuals.clear();
}

Individual *TemplatedPedigree::DrawIndividual(PoolManager &pools, TemplatedIndividual& tpl) {
	Individual *ind = NULL;
	IndID id= tpl.GetID();
	int statusFilter[] = {2,0,1,3};
	map<IndID, Individual*>::iterator membersEnd = members.end();
	if (members.find(id) != membersEnd)
		return members[id];
	if (tpl.IsFounder()) {
		ind = new Individual(tpl.id, tpl.pedID, pools.GetPoolCount());
		pools.DrawIndividual(*ind, tpl.gender == 2);
		ind->InitAlleleSource();
		ind->SetStatus(statusFilter[tpl.status]);
	}
	else {
		if (members.find(tpl.GetPatID()) != membersEnd && members.find(tpl.GetMatID()) != membersEnd) {
			Individual *father = members[tpl.GetPatID()];
			Individual *mother = members[tpl.GetMatID()];
			ind = father->DLCross(mother, tpl.GetID().indID, tpl.gender);
			members[id] = ind;
			//pools.DrawIndividual(*ind, tpl.gender == 2);
			ind->SetStatus(statusFilter[tpl.status]);
#ifdef DEBUG_VERBOSE
cerr<<"\n";
father->WritePhased(cerr, father->GetID(), false);
mother->WritePhased(cerr, mother->GetID(), false);
ind->WritePhased(cerr, ind->GetID(), false);
#endif

		}
		else
			cerr<<"Individual "<<tpl.pedID<<":"<<tpl.id<<" ("<<tpl.patID<<" - "<<tpl.matID<<") couldn't be rendered\n";
	}
	return ind;
}

//Technically, we should not have to do this, since the calling routine is responsible for the array of memory
void TemplatedPedigree::PurgeMembers() {
	members.clear();
}


void TemplatedPedigree::InitFounders(PoolManager& pools) {
	PurgeMembers();
	PurgeIndividuals(founders);

	IndividualLookup::iterator iTempl = individuals.begin();
	IndividualLookup::iterator iEnd = individuals.end();
	while (iTempl != iEnd) {
		if (iTempl->second->IsFounder()) 
			founders.push_back(DrawIndividual(pools, *(iTempl->second)));
		iTempl++;
	}
}

void TemplatedPedigree::ApplyPresentGenotypes(PoolManager &pools, Individual &person) {
	PoolManager::Iterator itr = pools.GetIterator();
	ChromPool *pool = itr.GetNext();

	while (pool) {
		pool->ApplyPresentGenotypes(&person, true);
		pool = itr.GetNext();
	}
		
}
TemplatedPedigree::StatusCounts TemplatedPedigree::BuildSample(PoolManager& pools, vector<Individual*>& individuals, PenetranceModel* model, bool failOnError) {
	StatusCounts curSet;
	StatusCounts sumCounts;
	float ambigRatio = statusCounts.AmbiguousRatio();
	int attempts = 0;
	bool doContinue = true;		// We'll repeatedly build a new set until we get the status' that we want
	while (doContinue && attempts < maxAttempts) {
		totalAttempts++;
		PurgeIndividuals(individuals);
		vector<Individual*>::iterator itr = founders.begin();
		vector<Individual*>::iterator end = founders.end();
		while (itr != end) {
			Individual *ind = *itr++;
			IndID id(ind->GetPedigreeID(), ind->GetID());
			if (members.find(id) == members.end()) {
				Individual *i = ind->Clone();
				members[id]=i;
			}
		}
	
		queue<IndID> ids = origOrder;
		map<IndID, Individual*>::iterator membersEnd = members.end();
		while (!ids.empty()) {
			IndID &id = ids.front();

			//Just go to next, if this is in the members already (founder)
			if (members.find(id) == membersEnd) {
				TemplatedIndividual *iTempl = this->individuals[id];
				Individual *ind = DrawIndividual(pools, *iTempl);
				if (ind == NULL)
					ids.push(id);
			}
			ids.pop();
		}
		//OK, this time we are writing them to the vector...in order
		ids = origOrder;
		curSet.Reset();
		while (!ids.empty()) {
			IndID &id = ids.front();
			assert(members[id] != NULL);
			Individual *ind = members[id];
			ApplyPresentGenotypes(pools, *ind);
			if (ind->GetStatus() != 3) {			// 3 indicates unknown, otherwise, it's known or ambiguous
				if (ambigRatio > Utility::Random::globalGenerator.drand()) {
					ind->SetStatus(2);
					curSet.Append(0);
				}
				else {
					ind->ApplyStatus(model);
					curSet.Append(ind->GetStatus()+1);
				}
			}
			else
				curSet.Append(3);
			ids.pop();
			individuals.push_back(ind);
		}
		assert(individuals.size() == members.size());
		PurgeMembers();
		doContinue = !statusCounts.Evaluate(curSet, tolerance);
		sumCounts = sumCounts+curSet;
		attempts+=doContinue;
/*		if (doContinue && attempts >= maxAttempts && failOnError) {
			sumCounts.Report(cerr, "Mean Simulated");
			cerr<<"genomeSIMLA was unable to produce a dataset that sufficiently matched the template pedigree within the specified:\n";
			cerr<<"-- # Attempts: "<<maxAttempts<<"\n";
			cerr<<"-- Threshold:  "<<tolerance<<"\n";
			exit(1);
		}
*/	}
	//if (attempts>0 && failOnError)
	//	sumCounts.Report(cerr, "Mean Simulated");
	curSet.attempts = sumCounts.attempts;
	return curSet;
}


uint TemplatedPedigree::ResetAttemptCount() {
	uint ac = totalAttempts;
	totalAttempts = 0;
	return ac;
}
TemplatedIndividual *TemplatedPedigree::GetMember(uint pedID, uint indID) {
	IndID id = IndID(pedID, indID);
	TemplatedIndividual *ind = NULL;
	if (individuals.find(id) != individuals.end())
		ind = individuals[id];
	return ind;
}

size_t TemplatedPedigree::GetMemberCount() {
	return origOrder.size();
}


}

}
