//
// C++ Implementation: pedigreereferencesample
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreereferencesample.h"

namespace Simulation {

bool PedigreeReferenceSample::DoWriteLOD = false;

float PedigreeReferenceSample::minMaxTolerance = 1.0;	
uint PedigreeReferenceSample::maxAttemptsAtReplication = 10000;

PedigreeReferenceSample::PedigreeReferenceSample(int replicates, float gtError, float phError, float missingData, const char *desc) :
		PedigreeSample(gtError, phError, missingData, desc), replicates(replicates), curModel(NULL), pools(NULL), minAffectedCount(0), maxAffectedCount(0),writeGenotypesForAll(false)
{
}


PedigreeReferenceSample::~PedigreeReferenceSample()
{
}	


/**
 * @brief Apply disease to the pedigree
 */
void PedigreeReferenceSample::BuildSample(PoolManager &pools, PenetranceModel* model) {
	curModel = model;
	this->pools=&pools;
	pedigrees.InitFounders(pools);

	pedigrees.BuildSample(pools,people, model, false);
}
	
void PedigreeReferenceSample::ApplyPhenocopyError(PoolManager *pools, PenetranceModel *model) {
	//We haven't really discussed this yet
}
	
string PedigreeReferenceSample::GetDescriptor() {
	stringstream ss;
	ss << "-ref" << description;
	ss<<".ped";
	return ss.str();
}


string PedigreeReferenceSample::GetSummary() {
	assert(0);
}

string PedigreeReferenceSample::GetDetails() {
	assert(0);
}

string PedigreeReferenceSample::Details() {
	return string("Pedigree Reference Sample: ") + description;
}

void PedigreeReferenceSample::SetReplicateCount(int replicateCount) {
	replicates = replicateCount;
}

PedigreeTemplates::TemplatedPedigree::StatusCounts PedigreeReferenceSample::Replicate() {
	Purge();
	PedigreeTemplates::TemplatedPedigree::StatusCounts counts = pedigrees.BuildSample(*pools, people, curModel);
	return counts;

}

int PedigreeReferenceSample::WriteDataset(const char *baseFilename, uint *gtCount) {
	string path, base, ext;
	SplitIntoComponents(baseFilename, path, base, ext);
	int observations = 0;
	if (path == "")
		path = ".";
	map<int, int> attempts;
	for (int i=0; i<replicates; i++) {
		stringstream ss;
		PedigreeTemplates::TemplatedPedigree::StatusCounts counts = Replicate();
		ss.str("");
		if (replicates > 1)
			ss<<path<<DIR_SLASH<<base<<"-r"<<i<<"."<<ext;
		else
			ss<<baseFilename;

		counts.Report(cerr, ss.str().c_str());
		int tries = counts.attempts;
		if (attempts.find(tries)==attempts.end())
			attempts[tries] = 1;
		else
			attempts[tries]++;

		if (counts.attempts < PedigreeTemplates::TemplatedPedigree::maxAttempts) {
			ofstream file(ss.str().c_str());
			observations =+ WriteDataset(file, gtCount);
			file.close();
		}
		else {
			cerr<<"Unable to make model fit requirements. No dataset written: "<<ss.str()<<"\n";
			exit(1);
		}
		ss.str("");
		if (replicates > 1)
			ss<<path<<DIR_SLASH<<base<<"-r"<<i<<".lod";
		else
			ss<<baseFilename<<".lod";

		if (DoWriteLOD)
			WriteLOD(ss.str().c_str());


	}
	uint attemptCount = pedigrees.ResetAttemptCount();
	cout<<replicates<<" replicates completed.\t"<<attemptCount-replicates<<" failures.\n";
	return observations;
}



/**
 * Evaluate R/NR based on the presence of XO between SNP and disease 
 * locus (R).
 */
void PedigreeReferenceSample::WriteLOD(const char *baseFilename) {
	ofstream file(baseFilename);
	
	vector<Individual*>::iterator pItr = people.begin();
	vector<Individual*>::iterator pEnd = people.end();
	int locusCount = (*pItr)->GetSnpCount(0);
	vector<double> distances(locusCount);
	int diseaseLocus = ((DiseaseModel*)curModel)->GetLocus(0).locusIdx;
	(*pItr)->ReferenceDistance(distances, diseaseLocus);
	vector<int> rCount(locusCount, 0);
	vector<int> nrCount(locusCount, 0);
	
file<<"Testing the LOD calculation. Here are some details that might help track down errors:\n";
file<<"Disease Locus: "<<diseaseLocus<<"\n";
file<<"Number of disease Locus: "<<curModel->GetModelSize()<<"\n";

	while (pItr!=pEnd){ 
//cerr<<"*"; cerr.flush();
		(*pItr)->CalculateLOD(rCount, nrCount, diseaseLocus);
		pItr++;
	}


for (size_t i=0; i<rCount.size(); i++)
	cerr<<i<<"\t"<<rCount[i]<<"\t"<<nrCount[i]<<"\n";


	cerr<<"    LOD           r\tt^r\t(1-t)^nr\tt\t: numerator/denom\n";
	double dl = distances[diseaseLocus];
	for (int i=0; i<locusCount; i++) {
		if (i<diseaseLocus)
			file<<i+1<<"\t"<<diseaseLocus+1<<"\t";
		else
			file<<diseaseLocus+1<<"\t"<<i+1<<"\t";
		double t = fabs(dl-distances[i]) / 100.0;
		if (t > 0.5)
			t = 0.5;
		if (t == 0.0)
			file<<"0.0\n";
		else {
			double r = rCount[i];
			double nr = nrCount[i];
			double num = (pow(t, r) * pow(1.0-t,nr));
			double denom = pow(0.5, r + nr);
			file<<setiosflags(ios::scientific | ios::showpoint)<<setprecision(8)<<log10(num/denom)<<"\t"<<r<<" ("<<pow(t, nr)<<")"<<"\t"<<r<<" ("<<pow(1.0-t, nr)<<") "<<"\t"<<t<<" : "<<num<<"/"<<denom<<"\n";
			cout<<setiosflags(ios::scientific | ios::showpoint)<<setprecision(8)<<log10(num/denom)<<"\t"<<r<<" ("<<pow(t, r)<<")"<<"\t"<<nr<<" ("<<pow(1.0-t, nr)<<") "<<"\t"<<t<<" : "<<num<<"/"<<denom<<"\n";
		}
	}


	//Fixed Thetas
	double theta[] = { 0.1, 0.2, 0.3, 0.4 };
	cerr<<"    Fixed Theta\n";
	cerr<<"    LOD           r\tt^r\t(1-t)^nr\tt\t: numerator/denom\n";
	for (int j=0; j<4; j++) {
		for (int i=0; i<locusCount; i++) {
			if (i<diseaseLocus)
				file<<i+1<<"\t"<<diseaseLocus+1<<"\t";
			else
				file<<diseaseLocus+1<<"\t"<<i+1<<"\t";
			
			double t = theta[j];
			{
				double r = rCount[i];
				double nr = nrCount[i];
				double num = (pow(t, r) * pow((1.0-t),nr));
				double denom = pow(0.5, r + nr);
				file<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(8)<<log10(num/denom)<<"\t"<<r<<" ("<<pow(t, nr)<<")"<<"\t"<<r<<" ("<<pow(1.0-t, nr)<<") "<<"\t"<<t<<" : "<<num<<"/"<<denom<<"\n";
				cout<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(8)<<log10(num/denom)<<"\t"<<r<<" ("<<pow(t, r)<<")"<<"\t"<<nr<<" ("<<pow(1.0-t, nr)<<") "<<"\t"<<t<<" : "<<num<<"/"<<denom<<"\n";
			}
		}
	}
}

void PedigreeReferenceSample::PrepSample(PoolManager *pools, PenetranceModel* model) {
	//Load reference
	int individualCount = pedigrees.Load(reference.c_str());
	if (individualCount < 1) {
		cerr<<"There were no individuals present after loading the pedigree template: "<<reference<<". Please make sure this file is correctly formatted as a pedigree file with 6 column headers. \n";
		exit(1);
	}
	//I think the founder initialization should occur elsewhere
	//pedigrees.InitFounders()
}

void PedigreeReferenceSample::SetReference(const char *filename) {
	reference = filename;
}


int PedigreeReferenceSample::WriteDataset(ostream &os, uint *gtCount) {
	//vector<FamilyPositionNode> startPoints = familyStartPositions;
	//random_shuffle(startPoints.begin(), startPoints.end(), Utility::Random::globalGenerator);
	pools->ResolveGenotypes(people, curModel);

	uint count = people.size();
	
	uint observations = 0;
	for (uint i=0; i<count; i++)  {
		Individual *person=people[i];

		if (person->DoIncludeInDataset()) {
			observations++;
			PedigreeTemplates::TemplatedIndividual *tpl = pedigrees.GetMember(person->GetPedigreeID(), person->GetID());
			if (tpl->hasGenotypes || writeGenotypesForAll)
				person->WritePedigree(os, gtCount);
			else
				person->WritePedigreeMetaData(os);
		}
	}
	return observations;
}
uint PedigreeReferenceSample::GetCount() {
	return pedigrees.GetMemberCount();
}


}
