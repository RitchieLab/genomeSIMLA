//
// C++ Implementation: modellist
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "modellist.h"
#include <fstream>
#include <iostream>
#include "diseasemodeldetails.h"
namespace GenomeSIM {

namespace GUI {
int 		ModelList::instanceCount 	= 0;
ModelList  *ModelList::instance			= NULL;


void ModelList::Load(const char *filename) {
	ifstream file(filename);

	this->filename = filename;
	
	if (file.is_open()) {
		isLoading=true;
		char curLine[4096];
		while (!file.eof()) {
			file.getline(curLine, 4096);
			stringstream line(curLine);
			string label = "", sourcefile="", modelType="";
			bool writeEnabled;
			int locusCount;
			line>>label>>modelType>>writeEnabled>>locusCount>>sourcefile;
			while (!line.eof()){
				string word;
				line>>word;
				sourcefile+=" "+word;
			}

			AddFile(label.c_str(), modelType.c_str(), sourcefile.c_str(), writeEnabled, locusCount);
		}
		isLoading=false;
	}
	else
		cout<<filename<<" was not opened\n";
}

void ModelList::Save() {

	cout<<"Trying to save model list: "<<filename<<"\n";
	ofstream file(filename.c_str());
	vector<ModelFileItem>::iterator itr = contents.begin();
	vector<ModelFileItem>::iterator end = contents.end();

	while (itr != end) {
		if (!itr->toBeDeleted)
			file<<itr->label<<"\t"<<itr->modelType<<"\t"<<itr->editable<<"\t"<<itr->locusCount<<"\t"<<itr->filename<<"\n";
		itr++;
	}
}

int ModelList::AddFile(const char *label, const char *modelType, const char *filename, bool editable, int locusCount) {
	if (strlen(label) > 0 && strlen(filename) > 0) {
		contents.push_back(ModelFileItem(label, modelType, filename, editable, locusCount));
	}
	if (!isLoading)
		Save();
	return contents.size();
}


void ModelList::MarkForDeletion(int idx) {
	assert((size_t)idx < contents.size());
	contents[idx].toBeDeleted = true;
}

ModelList::ModelFileItem &ModelList::GetFile(int idx) {
	assert((size_t)idx < contents.size());
	return contents[idx];
}

void ModelList::Update(int idx, ModelList::ModelFileItem& item) {
	assert(idx < contents.size());

	DiseaseModelDetails *model = DiseaseModelDetails::OpenDiseaseModel(item.filename.c_str());

	if (model) {
		model->Load();

		item.modelType = model->GetType();
		item.locusCount = model->GetLocusCount();

		delete model;
	}
	contents[idx] = item;

	Save();
}	

string ModelList::GetFilename(int idx) {
	stringstream val;
//	string val="";
	if (idx<contents.size())
		val<<contents[idx].filename;
//		val=contents[idx].filename;
	
	return val.str();
}

int ModelList::GetFileCount() {
	return contents.size();
}



}

}
