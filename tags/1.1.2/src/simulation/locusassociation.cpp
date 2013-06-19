//
// C++ Implementation: locusassociation
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locusassociation.h"
#include "utility/strings.h"
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Utility;

namespace Simulation {

LocusAssociation::LocusAssociation(int pos, int level) : position(pos), level(level), 
			probability((int)pow((float)2, (float)level), 0.0), observedZeros(0), totalObservations(0) { }

LocusAssociation::~LocusAssociation() { }

bool LocusAssociation::GetAllele(Utility::Random& r) {
	int mul = 1;
	int index = 0;
//cout<<"Precedent: "<<precedent<<"\n";
	for (int i=level-1; i>=0; i--) {
		index+=precedent[i]*mul;
		mul*=level;
	}
	bool isOne = r.drand() < probability[index];
	if (!isOne)
		observedZeros++;
	return isOne;
}

bool LocusAssociation::Load(istream& input) {
	string snpID;
	input>>snpID;
	int count = probability.size();
	for (int i=0; i<count; i++) {
		float &cur=probability[i];
		input>>cur;
	}
	return true;
}

float LocusAssociation::GetFrequencyZero() {
	return (float)observedZeros/(float)totalObservations;
}

void LocusAssociation::ResetObservations() {
	observedZeros = totalObservations = 0;
}


void LocusAssociationGrid::ResetObservations() {
	vector<LocusAssociation>::iterator itr = loci.begin();
	vector<LocusAssociation>::iterator end = loci.end();

	while (itr != end) {
		itr->ResetObservations();
		itr++;
	}
}

BitSetType LocusAssociationGrid::GenerateChromosome(Utility::Random& rnd) {
	BitSetType chrom(loci.size(), false);
	BitSetType precedent(gridDepth, false);
	
	for (int i=0; i<gridDepth; i++) {
		if (rnd.drand() > 0.5)
			precedent[i]=true;
	}
	bool currentValue;
	vector<LocusAssociation>::iterator itr = loci.begin();
	vector<LocusAssociation>::iterator end = loci.end();
	int idx = 0;
	while (itr != end) {
		itr->precedent = precedent;
		currentValue = itr->GetAllele(rnd);
		chrom[idx++]=currentValue;
		precedent<<=1;
		precedent[0]=currentValue;
		itr++;
	}
	return chrom;
}

int LocusAssociationGrid::Load(const char *filename) {
	ifstream file(filename);
	
	if (file.is_open()) {
		char line[65536];
		file.getline(line, 65536);			///<This is header junk
		gridDepth = CountColumns(line);	///<First column is header
//		if (gridDepth > 1)
			gridDepth = (int)log2(gridDepth);

		///Lets make sure that the number of columns are correct
		//assert(gridDepth * gridDepth == cols);

		cout<<"Loading Association Grid:\n";
		while (!file.eof()) {
			file.getline(line, 65536);
			if (line[0] != '#' && line[0]!= '\n') {
				stringstream ss(line);
				LocusAssociation assoc(loci.size(), gridDepth);
				if (assoc.Load(ss))
					loci.push_back(assoc);
			}
			cout<<"*";
		}
		cout<<"\n";
	}
	else {
		cout<<"Unable to open Locus Association Grid: "<<filename<<"\n";
		exit(1);
	}
	return loci.size();
}

}
