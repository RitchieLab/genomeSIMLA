//
// C++ Implementation: locusmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locusmanager.h"
#include <iostream>
#include <fstream>
#include <wx/filedlg.h>
#include <wx/filename.h>

namespace GenomeSIM {

namespace GUI {


void FileBasedChromosome::Purge() {
	vector<Locus *>::iterator itr = loci.begin();
	vector<Locus *>::iterator end = loci.end();

	while (itr != end) {
		delete *itr;
		itr++;
	}
	loci.clear();
}

void FileBasedChromosome::Save() {
	Save(filename.c_str());
}

FileBasedChromosome::~FileBasedChromosome( ){
	Purge();
}

void FileBasedChromosome::Save(const char *filename) {
	ofstream file(filename);
	file.setf(ios::fixed, ios::floatfield);
	
	int width=16;

	size_t lociCount = loci.size();

	file << "Locus Log for "<<label<<"\n";
	file << lociCount<<" Loci\n";

	if (lociCount > 0) 
		loci[0]->WriteHeader(file, width);
	
	for(uint i=0; i<lociCount; i++)
		loci[i]->WriteFormatted(file, width);
	
}



bool FileBasedChromosome::Load(const char *filename) {
	this->filename = filename;
	return Load();
}

bool FileBasedChromosome::Load() {
	ifstream file(filename.c_str());


	if (!file.is_open()) {
		cout<<"There was a problem opening the file: "<<filename<<"\n";
		status = Error;
		return false;
	}

	char line[4096];
	//Eat the first three rows
	file.getline(line, 4096);
	file.getline(line, 4096);
	file.getline(line, 4096);

	int locCount = 0;
	Purge();

	Locus *loc;
	while (!file.eof()) {
		loc = new Locus(chrID, locCount);
		file>>*loc;	
		if (loc->Freq1() + loc->Freq2() > 0.0) {
			loci.push_back(loc);
			locusLookup[loc->GetLabel()] = locCount;
			locCount++;
		}
		else
			delete loc;
	}
	status = Open;
	return true;
}

bool FileBasedChromosome::HasValidFile() {
	bool hasLoci = loci.size() > 0; 
	if (status == Unknown) 
		return hasLoci && filename != ""; 
	else if (status == Open || status == Closed) 
		return hasLoci; 
	else 
		return false; 
}


Locus *FileBasedChromosome::GetLocus(const char *lbl) { 
	cout<<locusLookup.size()<<"\n";
	if (locusLookup.find(lbl) != locusLookup.end())
		return loci[locusLookup[lbl]]; 
	else
		return NULL;
}

bool FileBasedChromosome::RequestFile() {
	bool success = false;
	static string lastChrom("");		///Just used to resume the last directory they were looking in
	wxFileDialog chromSelect(NULL, wxT("Select Chromosome"), wxT(""), wxT(lastChrom.c_str()), wxT("Locus Report(*.loc)|*.loc|All Files(*.*)|*.*|Affy Files(*.affy)|*.affy"), wxOPEN);
	if (chromSelect.ShowModal() == wxID_OK) {
		wxString path = chromSelect.GetPath();
		wxFileName relFilename = path;
		relFilename.MakeRelativeTo(wxGetCwd());
		lastChrom = filename = relFilename.GetFullPath().c_str();
		
		success = Load();
		
	}
	return success;
}

vector<Locus> FileBasedChromosome::ExtractLoci() {
	vector<Locus> loc;
	size_t locCount = loci.size();
	for (size_t i=0; i<locCount; i++) {
		loc.push_back(*(loci[i]));
	}
	return loc;
}

}
}
