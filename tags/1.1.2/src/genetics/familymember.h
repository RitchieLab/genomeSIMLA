//
// C++ Interface: familymember
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICSFAMILYMEMBER_H
#define GENETICSFAMILYMEMBER_H
#include "utility/utility.h"
#include <map>
#include "genotypedata.h"

namespace Genetics {
using namespace std;
using namespace Utility;



/**
@brief Represents a single entity from within a family. 
This will aid in the identification of various relationships required for cleaning data prior to certain pedigree evaluations

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class FamilyMember  {
public:
	

class GenoLkup : public TypeLookup {
public:
	GenoLkup(int badIdx=0) {
		this->badIdx = badIdx;
	};
	~GenoLkup() {};
	string GetValue(int i);
	int GetValue(const char* key);				///<Returns the mapped value setting it if necessary
	ostream *Report(ostream* os);					///<Appends the mapping to the stream
	const char *GetHeader();
	ostream &operator>>(ostream& os);	
		
}; 

	/**
	 * @brief Primary Constructor
	 * @param father Pointer to the family member IDd as father
	 * @param mother Pointer to the family member IDd as mother
 	 * @param genotypeCount Number of genotypes to be associated with the snp
	 * @param id Identification ID
	 */
    FamilyMember(FamilyMember *father, FamilyMember *mother, const char *famID, const char *id);//, uint snpID);
	FamilyMember(const char *famID, const char *id);//, uint snpID);
    ~FamilyMember();

	/**
	 * Adds the child to the list of children
	 */
	void AddChild(FamilyMember *child);

	/**
	 * @brief Return the ID of the member
	 */
	const char *GetID();

	/**
	 * @brief Returns the families ID in character format
	 */
	const char *GetFamilyID();

	/**
	 * @brief Set parents 
	 */
	void SetParents(FamilyMember *father, FamilyMember *mother);

	/**	
 	 * @brief Return the parent pointers
	 */
	FamilyMember *GetFather();
	/**	
 	 * @brief Return the parent pointers
	 */
	FamilyMember *GetMother();

	/**
	 * @brief Returns the number of children associated with this node
	 */
	uint ChildCount();

	/**
	 * @brief Use this to evaluate affected status in case the value for affected changes
	 */
	bool IsAffected();

	/**
	 * @brief Returns the genotypedata, instantiating it, if it isn't in memory already. 
	 * This currently uses a constant, fixed, gt parser. This should be opened up when used for real
	 */
	GenotypeData *GetGenotypeData();
	GenotypeData *GetGenotypeData(string sibID);
	/**
	 * @brief Returns the genotype label of at a given position
	 */
	string GetGenotypeLabel(uint position);

	/**
	 * @brief Returns the genotype value at a given position
	 */
	uint GetGenotypeValue(uint position);
	uint GetGenotypeValue(uint position, string sibID);

	/**
	 * @brief Returns Genotype description of the individual at position, idx. 
	 * This depends on genetic data being available locally. If we rely on genetic
	 * data to be kept solely in the snp repository, we will need to update it!
	 * @param idx The 0 based index of the genotype of interest
	 * @return Allelic data associated the genotype at position idx
	 */	 
	string GetGtValue(uint idx);

	/**
	 * @brief Determine if an individual has been marked for deletion
	 * @return True/False to indicate that the indiviudal should be deleted
	 */
	bool MarkForDeletion();

	/**
	 * @brief Mark an individual for deletion (or not)
	 * @param doDelete True/False for the deletion status
	 */
	void MarkForDeletion(bool doDelete);
	
	/**
	 * @brief Have an individual dump it's data to a stream
	 */
	ostream *Report(ostream *os);
	
	/**
	 * @brief Set the affected status for an individual
	 */
	void SetStatus(int affectedstatus);

	/**
	 * @brief Set the id for the member's position in the SNP bit vectors
	 */
	void SetIndividualIdx(uint val);

	/**
	 * @brief Returns the index (into the snp genotype vectors) associated with the local individual
	 */
	uint GetIndividualIdx(uint idx);

	/**
	 * @brief Used for reporting the status of the family member to logs
	 */
	string inclusionStatus;

	/**
	 * @brief Set's the weight of an entry so it will be replicated according to the statistic
	 */
	void SetWeight(uint weight);

	/**
	 * @brief Configure an masked version of the genotype data so that we can account for missing data
	 * for any of the immediate family's sibs. 
	 */
	void BuildMaskedGT(FamilyMember *other);

	/**
	 * @brief Return the desired weight
	 */
	uint GetWeight();

	/**
	 * @brief append this member's status to the set of statuses
	 */ 
	void AppendStatus(CaseControlStatus &status);

	/**
	 * @brief We don't want to lose members, unless they are supposed to be "children"- but, we need to know
	 * if their status is unknown	
	 */
	bool IsUnknownStatus();

	void SetEffectiveStatus(int newStat) {
		effectiveStatus = newStat;
	}
	int GetEffectiveStatus() {
		return effectiveStatus;
	}
	void ResetEffectiveStatus() {
		effectiveStatus = affectedStatus;
	}
	
	/**
	 * Some global information relevant to all members
	 */
	static int _affectedValue;					///<The value used to denote affected status
	static int _unaffectedValue;				///<The value used to denote unaffected status

protected:
	FamilyMember *father;						///<Pointer to the father
	FamilyMember *mother;						///<The pointer to the mother 
	
	string famID;								///<Family ID used to associate an extended family
	string id;									///<The individual's id (from ped file)

	vector<FamilyMember *> children;			///<children associated with this member
	int affectedStatus;							///<Stores the affected status
	int effectiveStatus;						///<Sets the "effective" status based on permuted statuses (if applied)
	GenotypeData *genotypedata;					///<Pointer to the genetic data associated with this member
	map<string, GenotypeData*> maskedGT;		///<The genotype data masked to correspond with a sibling's missing data
	static GenoLkup lkupTable;					///<Used for translating the genotype data 
	bool markForDeletion;						///<Used to mark for purging
	vector<uint> individualIdx;					///<Since this can be weighted, ind Index is a list of indices							
	uint weight;								///<This is used by the parser to duplicate weighted entries

};

inline
void FamilyMember::SetIndividualIdx(uint val) {
	individualIdx.push_back(val);
}

inline
uint FamilyMember::GetIndividualIdx(uint idx) {
	assert(0);
 	assert(idx < individualIdx.size());
	return individualIdx[idx];
}

inline
uint FamilyMember::GetGenotypeValue(uint position) {
	if (genotypedata == NULL)
		return lkupTable.GetNotEncodedIdx();
	else
		return genotypedata->GetGenotypeIndex( position );
}

inline
GenotypeData *FamilyMember::GetGenotypeData(string sibID) {
	return maskedGT[sibID];
}

inline
uint FamilyMember::GetGenotypeValue(uint position, string sibID) {
	GenotypeData *gt=maskedGT[sibID];
	assert(gt != NULL);
	return gt->GetGenotypeIndex(position);
}

inline
string FamilyMember::GetGenotypeLabel(uint position) {
	if (genotypedata == NULL)
		return lkupTable.GetValue(lkupTable.GetNotEncodedIdx());
	else
		return genotypedata->GetGenotypeLabel( position );
}

inline
string FamilyMember::GetGtValue( uint idx) {
	return lkupTable.GetValue(idx);
}	

inline
uint FamilyMember::GetWeight() {
	return weight;
}

inline
void FamilyMember::SetWeight(uint weight) {
	this->weight=weight;
}

inline
void FamilyMember::AppendStatus(CaseControlStatus &status) {
	uint indexCount=individualIdx.size();
	uint value=0;

	//Going backwards to avoid resizing each time
	for (uint i=indexCount; i>0; i--) {
		value=individualIdx[i-1];
		status.SetStatus(value, IsAffected());

	}
}


inline
bool FamilyMember::MarkForDeletion() {
	return markForDeletion;
}

inline
void FamilyMember::MarkForDeletion(bool doDelete) {
	markForDeletion=doDelete;
}

inline
GenotypeData *FamilyMember::GetGenotypeData() {
	if (genotypedata == NULL)
		genotypedata=new GenotypeData(&lkupTable);
	return genotypedata;
	
}


inline
void FamilyMember::BuildMaskedGT(FamilyMember *other) {
	assert(genotypedata != NULL);
	//assume that this hasn't been built before
	GenotypeData *gt = genotypedata->BuildMaskedGT(other->genotypedata);
	maskedGT[other->GetID()]=gt;
}
/**
	@brief Manages the Family members
 */
class FamilyMemberPool {
public:
	/**
	 * @brief Return an instance tot he pool
	 */
	static FamilyMemberPool *Instance();

	/**
	 * @brief Release the instance of the pool
	 */
	static void Release();
	
	/**
	 * @brief Returns a pointer to a valid entry based on the individual's ID
	 */
	FamilyMember *GetFamilyMember(const char *famID, const char *id);	
	
protected:
	static FamilyMemberPool *_instance;			///<Singleton instance
	static int _instanceCount;					///<Singleton instance count

	map<string, FamilyMember *> individuals;	///<Data storage

	FamilyMemberPool() {};						///<Only this class can instantiate itself
	~FamilyMemberPool();						///<Only this class can delete itself
};

inline
FamilyMemberPool::~FamilyMemberPool() {
	map<string, FamilyMember *>::iterator end=individuals.end();
	map<string, FamilyMember *>::iterator itr=individuals.begin();
	
	for (; itr!=end; itr++)	{
		delete itr->second;
	}

	individuals.clear();
}

inline
FamilyMemberPool *FamilyMemberPool::Instance() {
	if (_instance == NULL) {
		_instance=new FamilyMemberPool();
		_instanceCount = 0;
	}
	_instanceCount++;
	return _instance;
}

inline
void FamilyMemberPool::Release() {
	
	if (--_instanceCount < 1) {
		delete _instance;
		_instance=NULL;
	}
		
}	

inline
FamilyMember *FamilyMemberPool::GetFamilyMember(const char *famID, const char *id) {
	string indID=string(famID) + "x" + string(id);
	map<string, FamilyMember *>::iterator end=individuals.end();
	map<string, FamilyMember *>::iterator start=individuals.find(indID.c_str());
	FamilyMember *individual;

	//If we didn't find them, let's add one
	if (start==end) {
		individual = new FamilyMember(famID, id);	//, indCount);
		individuals[indID]=individual;
	} else
		individual=start->second;
	
	return individual;
}

inline
FamilyMember::FamilyMember(const char *famID, const char *id) : father(NULL), mother(NULL), famID(famID), id(id), affectedStatus(0), effectiveStatus(0), genotypedata(NULL), markForDeletion(true) {}

inline
FamilyMember::FamilyMember(FamilyMember *father, FamilyMember *mother, const char *famID, const char *id) : father(father), mother(mother), famID(famID), id(id), affectedStatus(0), effectiveStatus(0), genotypedata(NULL), markForDeletion(true) {}

inline
FamilyMember::~FamilyMember() {
	if (genotypedata)
		delete genotypedata;
	map<string, GenotypeData*>::iterator end=maskedGT.end();
	map<string, GenotypeData*>::iterator itr=maskedGT.begin();
	
	for (; itr!=end; itr++) {
		delete itr->second;
	}

}

inline
void FamilyMember::SetParents(FamilyMember *father, FamilyMember *mother) {
	this->father=father;
	this->mother=mother;

	if (father)
		father->AddChild(this);
	if (mother)
		mother->AddChild(this);

}
inline
bool FamilyMember::IsUnknownStatus() {
	return (effectiveStatus != _affectedValue && effectiveStatus != _unaffectedValue);
	//return (affectedStatus != _affectedValue && affectedStatus != _unaffectedValue);
}
inline
bool FamilyMember::IsAffected() {
	return effectiveStatus == _affectedValue;
	//return affectedStatus==_affectedValue;
}


inline
void FamilyMember::SetStatus(int status) {
	if (status != _affectedValue && status != _unaffectedValue) {
		markForDeletion=true;
	} 
	effectiveStatus = affectedStatus=status;
}

inline
uint FamilyMember::ChildCount() {
	return children.size();
}

inline
FamilyMember *FamilyMember::GetFather() {
	return father;
}

inline
FamilyMember *FamilyMember::GetMother() {
	return mother;
}

inline
void FamilyMember::AddChild(FamilyMember *child) {
	children.push_back(child);
}

inline
const char *FamilyMember::GetID() {
	return id.c_str();
}

inline
const char *FamilyMember::GetFamilyID() {
	return famID.c_str();
}


inline
string FamilyMember::GenoLkup::GetValue(int i) {
	string value;
	if (i == 0)
		value="0 0";
	else if (i==1)
		value="1 1";
	else if (i==2)
		value="1 2";
	else if (i==3)
		value="2 2";
	else
		value="0 0 ";
	return value;
}
inline 
int FamilyMember::GenoLkup::GetValue(const char *key) {
	if (strncmp(key, "1 1",3)==0)
		return 1;
	else if (strncmp(key, "1 2",3)==0 || strncmp(key, "2 1",3)==0)
		return 2;
	else if (strncmp(key, "2 2",3)==0)
		return 3;
	else	{
		return 0;
	}
}

inline
ostream *FamilyMember::GenoLkup::Report(ostream *os) {
	*os<<"1 - 1 1\n2 - 1 2 or 2 1\n3 - 2 2\n0 - Anything else\n";
	return os;
}

inline
const char *FamilyMember::GenoLkup::GetHeader() {
	return headerName.c_str();
}

inline
ostream &FamilyMember::GenoLkup::operator>>(ostream &os) {
	Report(&os);
	return os;
}	



}

#endif
