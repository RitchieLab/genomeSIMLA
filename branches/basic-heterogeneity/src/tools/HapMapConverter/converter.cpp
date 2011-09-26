//
// C++ Implementation: genomesim
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "converter.h"
#include <iostream>
#include <iomanip>
#include "utility/lineparser.h"
#include "utility/strings.h"
#include "utility/executionlog.h"

#ifndef TEST_APP
int main(int argc, char **argv) {
	Tools::Converter app;
	
	if (app.ParseCmdLine(argc, argv))
		app.Start();
	return 0;
}
	
#endif

namespace Tools {

using namespace std;

Converter::Converter()  {
	appname 	= "HM-Converter";
	appfunction = "HapMap to %TODO-NAME_FORMAT%";
	major	 	= 1;
	minor 	 	= 0;
	bugFixes 	= 0;
}


Converter::~Converter() {
	PurgeLoci();
	PurgeIndividuals();
}

void Converter::PurgeIndividuals() {
	vector<Individual *>::iterator itr = individuals.begin();
	vector<Individual *>::iterator end = individuals.end();

	while (itr != end) {
		delete *itr++;
	}
}

void Converter::PurgeLoci() {
	vector<LocusConverter *>::iterator itr = loci.begin();
	vector<LocusConverter *>::iterator end = loci.end();
	
	while (itr != end) {
		delete *itr++;
	}

	loci.clear();
}

/**
 * @brief Prints the help contents
 */
void Converter::PrintHelp() {
	PrintBanner();
	
	cout<<"\n--------------------------------------\n";
	
	configuration.GenerateReport(cout);

}

void Converter::ParseInclusionsList() {
	ifstream file(configuration.inclusionsList.c_str());
	int inclusionCount=0;
	string label;
	map<string, LocusConverter*>::iterator end = locusLookup.end();
	while (!file.eof()) {
		label="";
		file>>label;
		if (label.length() > 0){ 
			if (locusLookup.find(label) != end) {
				locusLookup[label]->doWriteToFile = true;
				inclusionCount++;
			}
		}
	}

	cout<<"Total Inclusions: "<<inclusionCount<<"\n";
}


void Converter::LoadSampleDetails() {
	ifstream file(configuration.sampleInputFile.c_str());
	int individualCount = individuals.size();
	
	string id;
	string gender;

	for (int i=0; i<individualCount; i++) {
		file>>id>>gender;
		individuals[i]->SetLabel(id.c_str());
		individuals[i]->SetGender(atoi(gender.c_str()));
	}	
	cout<<individuals[0]->GetLabel()<<" to "<<individuals[individualCount-1]->GetLabel()<<"\n";
}

void Converter::LoadHaplotypes() {
	//Load the snp data
	ifstream file(configuration.phaseInputFile.c_str());
	Individual *ind = new Individual(&loci);

	while (ind->Parse(file)) {
		individuals.push_back(ind);
		ind = new Individual(&loci);
	}
	delete ind;
	
	cout<<"\nPhased file parsed: "<<individuals.size()<<"\n";
	file.close();

}

void Converter::LoadLegend() {
	ifstream file(configuration.legendInputFile.c_str());
	
	PurgeLoci();

	char line[4096];
	//Get rid of the header
	file.getline(line, 4096);
	
	LocusConverter *locus = new LocusConverter();
	while (locus->Parse(file)) {
		loci.push_back(locus);
		locusLookup[locus->label]=locus;
		locus = new LocusConverter();
	}
	delete locus;

	cout<<"Legend Parsed: "<<loci.size()<<"\n";
	cout<<"Last SNP: "<<loci[loci.size()-1]->label<<"\n";
	file.close();

	
	
	
}



/**	
 * @brief Starts execution
 */
void Converter::Start() {
	configuration.GenerateReport(cout);

	//Load the Legend first
	LoadLegend();
	LoadHaplotypes();
	LoadSampleDetails();
	ParseInclusionsList();

	WriteData();
}

void Converter::WriteData() {
	cout<<"Writing Data ("<<configuration.destinationFilename<<")\n";
	ofstream file(configuration.destinationFilename.c_str());

	int lociCount = loci.size();
	int indCount  = individuals.size();
	file<<indCount<<"\n"<<lociCount<<"\nIND";
	for (int i=0; i<indCount; i++) 
		file<<", "<<individuals[i]->GetLabel();
	file<<"\n";	
	
	for (int l=0; l<lociCount; l++){ 
		LocusConverter *locus=loci[l];
		if (locus->doWriteToFile) {
			file<<locus->label;
			for (int i=0; i<indCount; i++) 
				file<<", "<<individuals[i]->GetGenotype(l);
			file<<"\n";
		}
	}
}

/**
 * Grab the pieces from the command line
 */
bool Converter::ParseCmdLine(int argc, char **argv) {
	bool success=true;

	if (argc != 2) {
		PrintHelp();
		return false;
	}
		
	configuration.SetConfig(argv[1]);
	Utility::LineParser lp('#');
	if (lp.Parse(argv[1], &configuration) > 0) {

	}
	else {
		cout<<"No valid parameters was found in file, "<<configuration.configurationFilename<<"\n";
		success=false;
	}
		
	return success;
}


}
