//
// C++ Interface: individual
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATIONINDIVIDUAL_H
#define SIMULATIONINDIVIDUAL_H

#include <string>
#include <vector>
#include <map>
#include "chromosome.h"
#include "utility/utility.h"
#include "diseasemodel.h"
#include "utility/strings.h"
#ifdef CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#endif
#ifdef USE_XY
#include "cpair.h"
#endif

namespace Simulation {

using namespace std;
using namespace Utility;

/**
base class for a person

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class Individual{
public:
	/**
	 * @brief Construction for arrays: Usage requires a call to Init() afterward to get everything set up properly
	 */
	Individual();
    Individual(uint id, uint pedID, uint chrCount=1);
    virtual ~Individual();

	/**
	 * @brief Initialize the basic components (if using default constructor)
	 */
	void Init(uint id, uint pedID, uint chrCount=1);
	void InitAlleleSource();

	uint GetID() const { return id; }

	uint GetPedigreeID() { return pedID; }

	/**
	 * Return the pool ID associated with the maternal and paternal Chromosomes. 
	 * This value is used to help avoid duplicate individuals being drawn.
	 */
	std::string GetPoolID(uint chromeID);

	/**
	 * Set an ID associated with a given chromosome
	 */
	void SetPoolID(uint chromeID, const char *id);

	/**	
	 * Set the chromosomal data
	 */
	void SetChromosomalData(uint chromeID, Chromosome &p, Chromosome &m);
	void SetChromosomalData(Chromosome *p, Chromosome *m);
	bool NeedsCompleteGenotypeInformation();
	void SetChromIDs(uint chromID, uint chID1, uint chID2);
	uint ChrID1(uint chromID) { return chrID1[chromID]; }
	uint ChrID2(uint chromID) { return chrID2[chromID]; }
	void ExpandModelLoci(uint chID, vector<uint>& modelLoci);
	/**
	 * @brief Pass a model to the individual to allow it to determine the status
	 * @param model This is the model that is to be used in assigning status
	 * @return T/F indicating affected status
	 */
	virtual int ApplyStatus(StatusModel::DiseaseModel* model);
	vector<uint> GetDiseaseGenotypes(StatusModel::DiseaseModel *model);
	/**
	 * @brief Returns the genotype for a given locus
	 */	
	virtual int GetGenotype(uint chromID, uint locus);

	/**
	 * @brief Returns true if individual is missing data at locus
	 */
	virtual bool MissingGenotype(uint chromID, uint locus);


	virtual int GetChromosomeCount() { return chromosomeCount; }

	virtual int GetSnpCount();

	float MinorAlleleFreq(uint chID, uint locID);

	/**
	 * @brief flip one of the alleles according the loci direction
	 */
	virtual bool ChangeGenotype(uint chrID, uint locus, int errorDir);
	virtual bool ClearLocus(uint chrID, uint loc);
	/**
	 * @brief Setup the affected value
	 */
	virtual void SetStatus(int isAffected);
	virtual void SetGender(int gender);
	virtual float GetOutcome();
	virtual void SetOutcome(float& outcome);

	void ReferenceDistance(vector<double>& distances, int diseaseLocus);
	void CalculateLOD(vector<int>& r, vector<int>& nr, int diseaseLocus);
	/**
	 * @brief Used to determine if an individual is affected
	 */
	virtual bool IsAffected();
	
	uint GetStatus();
#ifdef USE_XY
	virtual int GetGenotypeXY(uint locus);
	virtual bool ChangeGenotypeXY(uint locus, uint errorDir);
	virtual bool ClearLocusXY(uint loc);
//	virtual bool SetChromosomalData(uint chromID, Chromosome &p, Chromosome &m);
//	virtual bool MissingGenotype(uint chrID, uint locus);
	AlleleSource<LocusXY>* DLCrossXY();
	AlleleSource<LocusXY>* CrossXY();

	void SetXX(AlleleSource<LocusXY>* x, AlleleSource<LocusXY>* y, PAR_Region<LocusXY>* par);
	void SetXY(AlleleSource<LocusXY>* x, AlleleSource<LocusXY>* y, PAR_Region<LocusXY>* par);
	bool VerifyChromosomeIs(AlleleSource<LocusXY>* chromosome, bool isX);
	bool VerifyChromosomeIs(bool isX);
#endif

	void WriteChromosome(std::ostream& os, uint chrID);

	/**
	 * @brief this only matters if you want valid pedigrees to be written
	 */
	void SetPedigreeMeta(Individual *dad, Individual *mom);
	/**
	 * @Brief Cross over to create a new individual
	 * @param father the male
	 * @param id the new individual's id
	 * @param gender if this is zero, gender is random (if XY), otherwise, this should eventually be used to for the XO starting point for XY crossover
	 */
	Individual* DLCross(Individual *father, size_t id, int gender = 0);
	void SetDLPedigreeMeta(Individual *dad, Individual *mom);

	Individual *GetMother() { return mom; }
	Individual *GetFather() { return dad; }

	//This is for amish pedigrees where we have individuals who we don't want to write genotyping information for, but we want them to be present in the dataset
	void WritePedigreeMetaData(ostream&);
	void WritePedigree(ostream& os, uint *genotypeCounts, bool writeOutcome = false);
	void WritePedigree(std::ostream& meta, std::ostream&genotypes, bool writeMissingData, bool doWriteStatus =true);
	void WriteContinuous(ostream& os, uint *genotypeCounts, bool doWriteStatus = true);
	void WriteContinuous(std::ostream& meta, std::ostream&genotypes, bool  writeMissingData, bool doWriteStatus /*=true*/);
	void WriteMDR(ostream& os, uint *genotypeCounts, bool doWriteStatus = true);
	void WriteMDR(std::ostream& meta, std::ostream&genotypes, bool  writeMissingData, bool doWriteStatus /*=true*/);
	void ReadBinaryMdr(std::ifstream *genotypes, bool readMissingData, bool doWriteStatus, vector<uint> &chrSize);
	/**
	 * @brief Write phased data to stream
	 * @param os Stream to be written to
	 * @param indID individual's ID
	 * @param useFreqThresh Indicate if we should skip loci whose minor allele freq is
	 * 						below the threshold
	 */
	void WritePhased(ostream& os, uint indID, bool useFreqThresh);
	void WriteForMendel(std::ostream& meta, std::ostream& genotypes);
	/**
	 * @brief Produce a brand new individual with genetic information from both the mother and the father (the local is considered mommy)
	 */
	Individual *Cross(Individual *other, size_t id);
	
	/**
	 * @brief Produce a new chromosome array based on a crossing of the local chromosomal data
	 */
	Chromosome *Cross();
	Chromosome *DLCross();

	bool DoIncludeInDataset();
	void DoIncludeInDataset(bool doInclude);

	/**
	 * @brief Ordinarily, we want our data with 10 header columns. However, some apps require 6. False sets it to 6.
	 */
	static bool StandardPedigreeHeader;
	static bool PhasedPedigrees;
	static float FtoM_BirthRatio;
	/**
	 * @brief Threshold for datasets (minor allele frequency). Any that fail to meet this are dropped from the datasets
	 */
	static float minAlFreqThreshold;

	void ResolveGenotypes(vector<uint> &modelLoci);
	void ResolveGenotypes(uint chID, vector<uint> &modelLoci);

	/**
	 * @brief Return the number of chromosomes in the individual
	 */
//	uint GetChromosomeCount();
	
	/**
	 * @brief Returns the number of snps at a given chromosome
	 */
	size_t GetSnpCount(uint chID) { return matChrom[chID].LociCount(); }

	Individual *Clone();

	void PopulateChromosomePool(uint chID, std::vector<Chromosome>& pool);


	void WriteXOReport(std::ostream& os) const {
		//Just for debugging purposes:
		string patID = "x ";
		string matID = "x ";
		if (dad)
			patID = ToString((int)dad->GetID());
		if (mom)
			matID  = ToString((int)mom->GetID());
		for (uint chID=0; chID<chromosomeCount; chID++) {
			os<<pedID<<" "<<id<<"p "<<patID<<"\t";
//			patChrom[chID].WritePedFormat(os, 0.0, 0, patChrom[chID].LociCount());
			patChrom[chID].WriteXOPoints(os);
			os<<"\n";
			os<<pedID<<" "<<id<<"m "<<matID<<"\t";
//			matChrom[chID].WritePedFormat(os, 0.0, 0, matChrom[chID].LociCount());
			matChrom[chID].WriteXOPoints(os);
			os<<"\n";
		}
	}

	float EvaluateKinship(const Individual& other);


protected:

	/**
	 * Pedigree information
	 */
	uint id;								///<The individual id
	uint pedID;							///<pedigree id
	uint patID;							///<Paternal id
	uint matID;							///<Maternal id
	bool includeInDataset;			///<Allow the system to mark individuals to not be written to datasets

   unsigned int status;				///<Overall status to be reported in the file
	int gender;							///<Gender, if it has been assigned
	float outcome;						///<This is the raw status value, from continuous outcome (if appropriate)
	Chromosome *matChrom;				///<Mother's contribution
	Chromosome *patChrom;				///<Father's contribution
	uint chromosomeCount;				///<the number of chromosomes associated with the individual
	string *genePoolIDs;				///<Used so we can return an individual to the pool to be redrawn later if it couldn't be used for some reason

	Utility::BitSetType *missingData;

	Individual *mom;					///<Just so we can find our way back to the parents
	Individual *dad;					///<Just so we can find our way back to the parents

	std::vector<uint> chrID1;			///<Used for remote chromosome pool resolution
	std::vector<uint> chrID2;			///< ^ ^ ^ ^ ^ 
#ifdef USE_XY
	CPairXY chromXY;					///<Is responsible for compiling genotypes
	PAR_Region<LocusXY> *par;			///<Used to generate XO events in Males
#endif


};



#ifdef CPPUNIT
class IndividualTest : public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( IndividualTest );
	CPPUNIT_TEST( TestKinship );
	CPPUNIT_TEST_SUITE_END();
public:
    IndividualTest();
    ~IndividualTest();


	void setUp();
	void tearDown();
	
	void TestKinship();
protected:
	LocusArray loci;
	vector<Individual*> individuals;
	vector<Chromosome*> chromosomes;
};
#endif


}

#endif
