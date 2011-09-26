//
// C++ Implementation: genegenemodelarchive.h
//
// Description: Serializes gene gene models from binary file
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genegenemodelarchive.h"
#include <stdlib.h>
using namespace std;

namespace Biofilter {

GeneGeneModelArchive::iterator GeneGeneModelArchive::begin() {
	return iterator(source);
}

GeneGeneModelArchive::iterator GeneGeneModelArchive::end() {
	return iterator(source, true);
}

GeneGeneModelArchive::GeneGeneModelArchive(const char *tempFilename, uint modelInit, uint modelMax, bool binaryStorage)
	: models(tempFilename, modelInit, modelMax), binaryStorage(binaryStorage) { }

GeneGeneModelArchive::~GeneGeneModelArchive() { }

void GeneGeneModelArchive::Insert(GeneGeneModel& model) {
	assert(model.ImplicationIndex() > 0);
	//model.Write(cout, false);
	models.Insert(model);
}

void GeneGeneModelArchive::DisplayContents(ostream& os) {
	char *modelArchive=tempnam("/tmp", "models.archiveXXXXXX");

	models.Archive(modelArchive);
	WriteGeneArchive(os, binaryStorage);


	set<GeneGeneModel *, ggmScoreGT> modelStorage;


	fstream file(modelArchive, ios::in);
	char *singleArch=tempnam("/tmp", "singles.archiveXXXXXX");
	fstream singular(singleArch, ios::out|ios::binary);
	GeneGeneModel *model;

	while (file.good() && !file.eof()) {
		model = new GeneGeneModel();
		if (model->Load(file, binaryStorage)) {
			if (model->ImplicationIndex() > 1)
				modelStorage.insert(model);
			else
				model->Write(singular, binaryStorage);
		}
		else {
			delete model;
		}
	}
	file.close();
	singular.close();
	set<GeneGeneModel *, ggmScoreGT>::iterator itr = modelStorage.begin();
	set<GeneGeneModel*, ggmScoreGT>::iterator end = modelStorage.end();

	while (itr != end) {
		(*itr)->Write(os, false);
		delete *itr++;
	}
	modelStorage.clear();
	singular.open(singleArch, ios::in|ios::binary);
	while (singular.good() && !singular.eof()){
		GeneGeneModel *model = new GeneGeneModel();
		model->Load(singular, binaryStorage);
		model->Write(os, false);
		delete model;
	}
	unlink(modelArchive);
	unlink(singleArch);
}

map<float, uint> GeneGeneModelArchive::Archive(const char *geneArchive, const char *modelArchive) {
	models.Archive(modelArchive);
	WriteGeneArchive(geneArchive);


	set<GeneGeneModel *, ggmScoreGT> modelStorage;
	map<float, uint> modelCounts;

	//When reading the archived version (from the buffer), we have to assume binary
	ios_base::openmode mode = (ios_base::openmode)0;
	if (binaryStorage)
		mode = ios_base::binary;
	char *singleArch=tempnam ("/tmp", "singles");

	fstream file(modelArchive, ios::in | ios::binary);
	fstream singular(singleArch, ios::out|mode);
	GeneGeneModel *model;

	uint modelCount = 0;
	while (file.good() && !file.eof()) {
		model = new GeneGeneModel();
		if (model->Load(file, true)) {
			modelCount++;
			float ii = model->ImplicationIndex();
			if (modelCounts.find(ii) == modelCounts.end())
				modelCounts[ii]=model->EstimateModelCount();
			else
				modelCounts[ii]+=model->EstimateModelCount();


			if (ii > 1)
				modelStorage.insert(model);
			else
				model->Write(singular, binaryStorage);
		}
		else {
			delete model;
		}
	}
	file.close();
	file.open(modelArchive, ios_base::out | mode );
	if (binaryStorage) {
		uint zero=0;
		file.write((char*)&zero, sizeof(uint));
		file.write((char*)&modelCount, sizeof(uint));
	}
	else
		file<<modelCount<<"\n";

	singular.close();
	set<GeneGeneModel *, ggmScoreGT>::iterator itr = modelStorage.begin();
	set<GeneGeneModel*, ggmScoreGT>::iterator end = modelStorage.end();

	while (itr != end) {
		(*itr)->Write(file, binaryStorage);
		delete *itr++;
	}
	singular.open(singleArch, ios::in|mode);
	while (singular.good() && !singular.eof()){
		GeneGeneModel *model = new GeneGeneModel();
		if (model->Load(singular, binaryStorage))
			model->Write(file, binaryStorage);
		delete model;
	}

	unlink(singleArch);
	return modelCounts;
}

void GeneGeneModelArchive::WriteGeneArchive(ostream& os, bool useBinary) {
	RegionLookup::iterator itr = GeneGeneModel::geneLookup.begin();
	RegionLookup::iterator end = GeneGeneModel::geneLookup.end();

	while (itr != end) 
		(itr++)->second->Write(os, useBinary);

}

void GeneGeneModelArchive::WriteGeneArchive(const char *geneArchive) {
	ios_base::openmode mode = ios_base::out;
	if (binaryStorage)
		mode = ofstream::out|ofstream::binary;
	ofstream file(geneArchive, mode);
	WriteGeneArchive(file, false);
}



uint GeneGeneModelArchive::SummarizeModelCounts(map<uint, uint>& scores, map<uint, Region*>& regions) {
	string tempFilename = tempnam("/tmp", "models.gg");
	models.Archive(tempFilename.c_str());
	ifstream file(tempFilename.c_str(), ios::binary);
	uint modelCount = 0;
	GeneGeneModel curGG;
	while (file.good() && !file.eof()) {
		if (curGG.LoadBinary(file)) {
			uint localTotal = curGG.EstimateModelCount();
			modelCount+=localTotal;
			curGG.Write(cerr);
			scores[(uint)curGG.ImplicationIndex()]+=localTotal;
		}
	}
	unlink(tempFilename.c_str());
	return modelCount;
}

}
