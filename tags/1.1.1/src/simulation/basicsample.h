//
// C++ Interface: population
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONPOPULATION_H
#define SIMULATIONPOPULATION_H
#include "utility/types.h"
#include "chrompool.h"
#include "penetrancemodel.h"

namespace Simulation {
using namespace StatusModel;

class PoolManager;

/**
 * @brief Base class for all sample types
 * Each Sample object represents a single group of datasets. They produce the datasets according to their configuration
 * and are capable of writing data in the correct format
 * Samples draw individuals from the pool manager- so more than 1 sample can be drawn from a set of chromosome pools
 * Status and is applied through the application of a model manager (which ultimately obtains a model). 
 */
class Sample {
public:
	Sample(float genoError, float phenoError, float missingData) : genoError(genoError), phenoError(phenoError), missingData(missingData), overallAffectedCounts(NULL), overallUnaffectedCounts(NULL) { }

	virtual ~Sample() {	Purge(); }
	/**
	 * @brief Constructs a population drawn from the local pool
	 * @param pools Pool manager containing each of the pools
	 * @param model Disease model to be applied to the sample
	 */
	virtual void BuildSample(PoolManager &pools, PenetranceModel* model)=0;

	/**
	 * @brief Adds individuals to the local sample based on the model passed
	 * @param pools Pool manager containing each of the pools
	 * @param model Disease model to be applied to the sample
	 * @param makeup The unit containing affected/unaffected counts
	 */
	//virtual void BuildSample(PoolManager &pools, PenetranceModel* model, UnitMakeup& makeup)=0;		

	/**
	 * @brief Apply genotype error based on each locus setting
	 */
	virtual void ApplyGenocopyError(PoolManager *mgr);
	
	/**
	 * @brief Apply phenocopy error 
	 */
	virtual void ApplyPhenocopyError(PoolManager *mgr, PenetranceModel* model)=0;

	/**
	 * @brief Apply a percentage to the dataset data and erase genotypes at the positions
	 */
	virtual void ApplyMissingData(PoolManager *mgr);

	/**
	 * @brief Works through the various types of errors and applies them each on the dataset
	 * according to their predefined settings
	 */
	virtual void ApplyErrors(PoolManager *mgr);
	/**
	 * @brief Resets the local data back to a starting point (for starting a new sample)
	 */
	virtual void Reset()=0;

	virtual void Purge();
	
	void WriteBinaryHeader(ostream &genotypes, uint fileVersion, const char *type);
	virtual string Details() = 0 ;
	/**
	 * @brief Returns a simple segment used for generating file names. This will be something like: 15-85.cc for 15% affected case/control. The application is still required to generate the rest of the filename, however
	 */
	virtual string GetDescriptor()=0;

	virtual string GetType()=0;

	virtual string GetLabel()=0;
	
	virtual string GetDetails()=0;

	virtual string GetSummary()=0;

	Individual *GetIndividual(size_t idx) { assert(idx < people.size()); return people[idx]; }
	size_t GetIndividualCount() { return people.size(); }

	/**
	 * @brief Dump the contents of the sample population to the stream
	 * @param os The stream to which the sample will be written
	 */
	virtual int WriteDataset(ostream &os, uint *gtCounts)=0;
	virtual int WriteDataset(const char *filename, uint *gtCounts);
	virtual int WriteBinaryDataset(ostream &meta, ostream &genotypes) = 0;
//	virtual int LoadBinaryDataset(ostream &metadata, ostream &genotypes) = 0;
	/**
	 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
	 * @param os The stream to which the sample will be written	
 	 */
	//virtual void WritePhased(ostream &os)=0;

//	virtual void GenerateLDMap(ostream &os)=0;
	
	virtual void GenerateReport(ostream &os, uint padding) = 0;

	/**
	 * @brief This allows us to determine the presence of meta data based on sample type (binary only). 
	 */
	virtual bool UsesMeta() { return false; }

	/**
	 * @brief returns t/f regarding the compatibility between the model and the dataset
	 */
	virtual bool Verify(PenetranceModel* model) { if (model->IsContinuous()) { cerr<<"Unable to mix discrete datasets with continuous disease models.\n"; return false;} return true; }

	/**
	 * @brief Initiates the valid type for the dataset if descibed in the file's contents
	
	static Sample *LoadSample(ifstream &metaData, ifstream &genotypes);
	Sample *LoadBinarySample(ifstream &metaData, ifstream &genotypes) = 0;
	 */

	virtual int GetSnpCount();

	static Sample *ParseBinaryFile(const char *gtFilename, const char *metaFilename, uint fileVersion);
	virtual void LoadBinarySample(ifstream *genotypes, ifstream *meta, uint peopleCount, uint chromCount, vector<uint> &chrID) = 0;
	class Iterator {
	public:
		/**
		* @brief Returns a pointer to the current individual or null
		*/
		Individual *operator->();			

		/**
	 	 * @brief returns a reference to the current individual.
	 	 * @note This could cause problems if the iterator is no longer pointing to a valid object
		 */
		Individual &operator&();
		
		/**
	 	 * @brief advances the operator
		 */
		Iterator &operator++();
	
		~Iterator();									///<Destruction
		void Reset();									///<Reset the iterator
		/**
	 	 * @brief Basic contruction
	   	 */
		Iterator(vector<Individual *> *repo) : repository(repo), position(repo->begin()) { }			
	protected:	
		vector<Individual *> *repository;				///<The repository we are iterating through
		vector<Individual *>::iterator position;			///<The current position
		
	};

	virtual bool ResolveGenotypes(ChromPool *pool, size_t i, bool justLoadedChrom);
	virtual void AppendGenotypeCountsToReport(StatusModel::ModelLociArray loci);
	
	string BuildGenotypeLabel(uint *genotypes, uint modelSize);
	virtual void ReportGenotypeCounts(StatusModel::ModelLociArray loci, ostream& os);
	void BuildGenotypeLabels(vector<string>& labels, uint modelSize);
	string BuildGenotypeLabel(uint genotype, uint position);
	uint GetMultiplier(uint genotype, uint position, uint modelSize);
	uint BuildGenotypeIndex(vector<uint>& genotypes);
	void WriteConfiguration(ostream& file);

	virtual void WriteSnpDetails(const char *, PoolManager *pools)=0;
	/**
	 * @brief Draws individual from the pool	 
	 * @param id individual ID
	 * @param pediID pedigree id
	 * @param pools The source of the chromosomes
	 * @param model The disease model (in order to assign status)
	 * @param status if not -1, this will redraw (the individual) until we get the desired status
	 * @param gender 0 - male, 1 - female
	 */
	Individual *DrawIndividual(int id, int pedID, PoolManager *pools, PenetranceModel* model, int status, bool gender);
	virtual void PrepSample(PoolManager *pools, PenetranceModel* model) { }

	void PopulatePool(std::vector<Chromosome>& pool, uint chrID);

protected:
	vector<Individual*> people;				///<The actual population
	float genoError;						///<Genotype error associated with the local sample
	float phenoError;						///<Phenotype error
	float missingData;						///<Missing data percentage
	uint aggregateCount;					///<How many groups (families or case/control pairs) are accounted for

	uint *overallAffectedCounts;			///<Affected Counts used in reporting genotype distribution
	uint *overallUnaffectedCounts;			///<Unaffected counts ''        
	uint *overallOtherCounts;               ///

};


/**
@brief Oversees the production of case/control datasets and saves them to file. 
These classes depend upon a remote GenePool array object to extract their individuals from

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BasicSample : public Sample {
public:
	/**
	 * @brief Construction
	 * @param percentAffected Percentage of individuals in the dataset to be found affected
	 * @param individualCount The total number of individuals to be found in a single dataset
	 */
    BasicSample(uint affCount, uint unaffCount, uint midCount, float genoErr, float phenoErr, float missing, const char *desc);

    virtual ~BasicSample();

	/**
	 * @brief Constructs a population based on the genetic data in pools with status derived from a 
	 * model provided by models.
	 * @param pools Array of chromosome pools. Each pool contains genetic data to be used to create a single individual
	 * @param chromCount Size of the array above
	 * @param models is used to randomly assign status
	 * @param individualCount Size of the resulting collection of individuals
	 */
	virtual void BuildSample(PoolManager &pools, PenetranceModel* model);


	/**
	 * @brief Resets the local data back to a starting point (for starting a new sample)
	 */
	virtual void Reset();


	/**
	 * @brief Returns a simple segment used for generating file names. This will be something like: 15-85.cc for 15% affected case/control. The application is still required to generate the rest of the filename, however
	 */
	virtual string GetDescriptor();

	virtual string GetType() { return "C/C"; }

	string GetSummary() { char msg[1000]; sprintf(msg, "%d Individuals", (affCount + unaffCount + midCount)); return msg;}

	string GetLabel() { return description; }

	string Details();

	string GetDetails();


	void WriteSnpDetails(const char *, PoolManager *pools);

	/**
	 * @brief Dump the contents of the sample population to the stream
	 * @param os The stream to which the sample will be written
	 */
	virtual int WriteDataset(ostream &os, uint *gtCounts);
	virtual int WriteBinaryDataset(ostream &metadata, ostream &genotypes);
//	virtual int LoadBinaryDataset(ostream &metadata, ostream &genotypes);
	/**
	 * @brief Dump the contents of the sample population to the stream (in phased haploview format)
	 * @param os The stream to which the sample will be written	
 	 */
	//virtual void WritePhased(ostream &os);

	void GenerateReport(ostream &os, uint padding);

	void ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model);

	void LoadBinarySample(ifstream *genotypes, ifstream *meta, uint peopleCount, uint chromCount, vector<uint> &chrID);
	
protected:
	//float percentAffected;					///<Percentage of overall individuals who are affected
	uint affCount;
	uint unaffCount;
	uint midCount;
	string description;
//	uint individualCount;

};

/**
 * This is experimental, and has not been verified. Do not use this 
 */
class ContinuousSample : public BasicSample {
public:
	/**
	 * @brief Construction
	 * @param percentAffected Percentage of individuals in the dataset to be found affected
	 * @param individualCount The total number of individuals to be found in a single dataset
	 */
    ContinuousSample(uint leftTail, uint rightTail, uint midCount, float genoErr, float phenoErr, float missing, const char *desc);
	/**
	 * @brief Dump the contents of the sample population to the stream
	 * @param os The stream to which the sample will be written
	 */
	virtual int WriteDataset(ostream &os, uint *gtCounts);
	virtual int WriteBinaryDataset(ostream &metadata, ostream &genotypes);
	virtual void ReportGenotypeCounts(vector<StatusModel::DiseaseLocus> loci, ostream& output);
	void AppendGenotypeCountsToReport(StatusModel::ModelLociArray loci);

	void ApplyPhenocopyError(PoolManager *pools, PenetranceModel* model);
	bool Verify(PenetranceModel* model);

};

inline
BasicSample::BasicSample(uint affCount, uint unaffCount, uint midCount, float genotypeError, float phenoError, float missing, const char *desc) : 
		Sample(genotypeError, phenoError, missing), 
		affCount(affCount), unaffCount(unaffCount), midCount(midCount), description(desc) {
}

inline
ContinuousSample::ContinuousSample(uint leftTail, uint midCount, uint rightTail, float genoErr, float phenoErr, float missing, const char *desc) :
		BasicSample(leftTail, rightTail, midCount, genoErr, phenoErr, missing, desc) { 
	cout<<"Genotype Error:		"<<genoError<<"\nPhenocopy Error:		"<<phenoError<<"\nMissing Data:		"<<missingData<<"\n";
}

inline
BasicSample::~BasicSample() {
 }


inline
void BasicSample::Reset() {
	people.clear();
}

inline
int Sample::GetSnpCount() {
	int count = 0;
	
	if (people.size() > 0)
		count = people[0]->GetSnpCount();

	return count;
}

inline
string BasicSample::Details() {
	return string("Case/Control: ") + description;
}
inline
string BasicSample::GetDescriptor() {
	stringstream ss;
	ss.precision(0);
	if (genoError > 0.0)
		ss<<"-"<<setprecision(2)<<(genoError*100.0)<<"GE";
	if (phenoError > 0.0)
		ss<<"-"<<setprecision(2)<<(phenoError*100.0)<<"PE";
	if (missingData > 0.0)
		ss<<"-"<<setprecision(2)<<(missingData*100.0)<<"M";
	ss<<setprecision(2)<<"-"<<description<<".mdr";
	return ss.str();
}


inline
Sample::Iterator::~Iterator() { }	

inline
Individual *Sample::Iterator::operator->() {
	Individual *ptr = NULL;
	if (position != repository->end())
		ptr = *position;
	return ptr;
}
inline
Individual &Sample::Iterator::operator&() {
	assert(position != repository->end());
	return *(*position);
}
inline
Sample::Iterator &Sample::Iterator::operator++() {
	position++;
	return *this;
}
inline
void Sample::Iterator::Reset() { 
	position = repository->begin();
}


}

#endif
