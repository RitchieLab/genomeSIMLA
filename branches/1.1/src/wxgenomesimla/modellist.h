//
// C++ Interface: modellist
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIMODELLIST_H
#define GENOMESIM_GUIMODELLIST_H
#include <string>
#include <vector>
#include <assert.h>
#include "utility/exception.h"

namespace GenomeSIM {

namespace GUI {

using namespace std;

/**
	@brief Provides an interface to the global list of models. 
	@Note Based on a modified Singleton design pattern. In this case, we have 2 methods
 		 	for pointer aquisition. One must be called first, and will be responsible for
	 		setting the filename for the list. The next, is the one to be used. If either
			is called in the wrong way, NULL is returned
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class ModelList{
public:
	struct ModelFileItem {
		string label;				///<Label of the model	
		string filename;			///<Name of the model's source file
		string modelType;			///<PEN, SIMPEN, SIMLA
		bool editable;				///<1/0 indicating editable status
		int locusCount;
		bool toBeDeleted;			///<If true, the model will not be saved

		ModelFileItem() : 
			label(""), filename(""), modelType(""), editable(0), locusCount(-1),toBeDeleted(false) {}
		ModelFileItem(const char *label, const char *type, const char *filename, bool editable, int locusCount) : 
			label(label), filename(filename), modelType(type), editable(editable), locusCount(locusCount), toBeDeleted(false){ } 
	};


	/**
	 * @brief Construct the structure with the name of the file. This MUST be the first instance call
	 */
	static ModelList *Instance(const char *filename);

	/**
	 * @brief Return the one and only instance. This can only be called after a filename has been established
	 */
	static ModelList *Instance();

	/**
	 * @brief Returns the number of instances currently in use
	 */
	static int GetInstanceCount();

	/**
	 * @brief returns the number of files in the repository
	 */	
	int GetFileCount();

	ModelFileItem &GetFile(int idx);

	/**
	 * @brief mark an item for deletion. The item won't go away until we reload the file- but it will be marked as such and not written to the archive list
	 */
	void MarkForDeletion(int idx);

	/**
	 * @Brief Returns the filename at index, idx
	 */
	string GetFilename(int idx);

	/**
	 * @brief adds a file to the list
	 * @param label How the model is shown to the user
	 * @param type What kind: PENTABLE, SIMPEN, SIMLA
	 * @return total number of items contained
	 */
	int AddFile(const char *label, const char *type, const char *filename, bool editable, int locusCount);
	
	
	/**
	 * @brief Releases the current instance (will release memory, if its the last
	 */
	void Release();

	/**
	 * @brief Allow a client to replace an existing item with a new one
	 */
	void Update(int idx, ModelList::ModelFileItem& item);

	void Save();				///<Write the labels & files to the file (skipping those to be deleted)

protected:
	static int instanceCount;
	static ModelList *instance;
	string filename;
	bool isLoading;				///<Just to make sure we don't save while we are loading
	/**
	 * @brief loads the file, filename
	 * @note This is called by the initial instance() method
	 */
	void Load(const char *filename);

	//Singleton
    ModelList();
    ~ModelList();


	vector<ModelFileItem> contents;

};



inline
void ModelList::Release() {
	if (--instanceCount == 0) {
		delete instance;
		instance=NULL;
	}
}

inline
int ModelList::GetInstanceCount() {
	return instanceCount;
}

inline
ModelList *ModelList::Instance(const char *filename) {

	cout<<"Initializing Model List: "<<filename<<"\n";
	assert(instance == NULL);
	instance = new ModelList();
	instanceCount = 1;
	try {
		instance->Load(filename);
		cout<<"Model List successfully opened\n";
	//We don't really care to handle this. It will be overwritten once we add some models
	}	catch (Utility::Exception::FileNotFound& e) {
		cout<<"File not opened\n";
	}
	return instance;
}
inline
ModelList *ModelList::Instance() {
	assert(instance);
	instanceCount++;
	return instance;
}

inline
ModelList::ModelList() : filename(""), isLoading(false) { }

inline
ModelList::~ModelList() { }

}

}

#endif
