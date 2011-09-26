//
// C++ Interface: familyrepository
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_EVALUATIONFAMILYREPOSITORY_H
#define GENETICS_EVALUATIONFAMILYREPOSITORY_H
#include <map>
#include "familynode.h"

namespace Genetics {
namespace Evaluation {
class FamilyRepoEvaluation;
}
using namespace std;
using namespace Evaluation;

/**
 * @brief Storage for all families
 */
class FamilyRepository {
public:
/**
 * @brief This allows basic iteration over the contents of the repository. 
 * The iterator can become invalid if changes are made to the repository (deletions) 
 */ 
class Iterator {
friend class FamilyRepository;
public:
	/**
	 * @brief Returns a pointer to the current family and advances to the next
	 * @return pointer to a family or NULL if we are at the end
	 */
	FamilyNode *GetNext();
	~Iterator();
	void Reset();
protected:	
	Iterator(map <string, FamilyNode *> *repository);
	map <string, FamilyNode *> *repository;
	map<string, FamilyNode *>::iterator position;
	
};



	/**
	 * @brief [singleton] Acquire an instance. There is currently assumed to be one repository
	 */
	static FamilyRepository *Instance();
	
	 /**
	  * @brief [singleton] Acknowledge the that we are done with our pointer
	  */
	static void Release();

	/**
	 * @brief Acquire a family node
	 */
	FamilyNode *GetNode(const char *famID);

	 /**
	  * @brief Perform an evaluation over the whole repository
	  */
	void PerformEvaluation(FamilyRepoEvaluation *eval);

	/**
	 * @brief Returns the number of families inside the repository
	 */
	uint GetFamilyCount(); 

	/**
	 * @brief Returns the number of parental groups available from all pedigrees
	 */
	uint CountSibships();

	/**
	 * @brief Dumps all families and their members
	 */
	void Purge();

	/**
	 * @brief Generates an overview of the contents of the repository
	 */
	void GenerateReport(ostream *os);

	/**
	 * @brief Removes entries that have been marked as "deleted"
	 */
	void PurgeDeleted();

	/**
	 * @brief Acquire the iterator for the first family
	 * @note This iterator doesn't guarantee that it will traverse the families in the same order each time 
	 */
	Iterator GetIterator() {
		return Iterator(&entries);
	}

	/**
	 * @brief Set up the log where families will report details during load
	 */
	void SetPedLog(BasicLog *log);

protected:
	static FamilyRepository *_instance;					///<Singleton instance
	static int _instanceCount;							///<Singleton instance count
	FamilyRepository() : pedLog(NULL) {};				///<Default constructor blocked from normal use
	~FamilyRepository();
	map <string, FamilyNode *> entries;					///<The storage of family data
	BasicLog *pedLog;									///<Log used for reporting on pedigree details
public:


};

inline
void FamilyRepository::PurgeDeleted() {
	map <string, FamilyNode *>::iterator itr = entries.begin();
	map <string, FamilyNode *>::iterator end = entries.end();
	
	for (; itr!=end; itr++) {
		if (itr->second->MarkForDeletion() ){
			cout<<"Purging member: "<<itr->second->GetID()<<"\n";
			delete itr->second;
			entries.erase(itr);
		}
	}	
}


inline
void FamilyRepository::Iterator::Reset() {
	position = repository->begin();
}

inline
FamilyRepository::Iterator::Iterator(map <string, FamilyNode *> *repository) : repository(repository) {
	position = repository->begin();
}



inline
FamilyRepository::Iterator::~Iterator() { }

inline
FamilyNode *FamilyRepository::Iterator::GetNext() {
	FamilyNode *current = NULL;
	if (position != repository->end()) { 
		current=position->second;
		position++;
	}
	return current;
}

inline
uint FamilyRepository::CountSibships() {
	map <string, FamilyNode *>::iterator itr = entries.begin();
	map <string, FamilyNode *>::iterator end = entries.end();
	uint pgCount = 0;
	for (; itr!=end; itr++) {
		pgCount += itr->second->CountSibships();
	}	
	return pgCount;

}

inline
uint FamilyRepository::GetFamilyCount() {
	return entries.size();
}

inline
FamilyRepository::~FamilyRepository() {
	Purge();
}
	
inline
void FamilyRepository::Purge() {
	map <string, FamilyNode *>::iterator itr = entries.begin();
	map <string, FamilyNode *>::iterator end = entries.end();
	
	for (; itr!=end; itr++) {
		delete itr->second;
	}	
	entries.clear();
}

inline
void FamilyRepository::SetPedLog(BasicLog *log) {
	pedLog = log;
}

inline
FamilyNode *FamilyRepository::GetNode(const char *famID) {
	FamilyNode *newFamily;
	map <string, FamilyNode *>::iterator itr;
	itr=entries.find(famID);
	if (itr == entries.end()) {
		newFamily=new FamilyNode(famID, pedLog);
		entries[famID]=newFamily;
	}
	else {
		newFamily=itr->second;
	}
	return newFamily;
}


inline
FamilyRepository *FamilyRepository::Instance() {
	if (_instance == NULL) {
		_instance=new FamilyRepository();
		_instanceCount=0;
	}
	_instanceCount++;
	return _instance;
}

inline
void FamilyRepository::Release() {
	assert(_instanceCount>0);
	if (--_instanceCount < 1) {
		delete _instance;
		_instance = NULL;
	}
}

}


#endif
