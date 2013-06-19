//Model.h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////// 
//
// Class contains penetrance tables and 
// list of loci associated with the model.
// The list of loci is kept as the locus indexes
// as kept in the Chromosome static data.
//
/////////////////////////////////////////////////////////////////////


#ifndef __BLANKMODEL_H__
#define __BLANKMODEL_H__

#include <string>
#include <iostream>
#include <vector>
#include "diseasemodel.h"
#include "poolmanager.h"


/*******************smod file***********************/
#define LABEL_DISEASELOCI 	"DISEASELOCI"
#define LABEL_PENTABLE 		"PENTABLE"
#define LABEL_FREQ_THRESH 	"FREQ_THRESHOLD"
#define LABEL_FREQ 			"FREQ"


#include "diseasemodel.h"
#include "diseaselocus.h"

namespace Simulation {
namespace StatusModel {

using namespace std;


class PenetranceModel : public DiseaseModel {
public:
	PenetranceModel();

	~PenetranceModel();

	/**
	 * @brief Load the model at the local file, filename
	 */
	virtual void Load()=0;

	virtual void Refresh(PoolManager *pools)=0;

	/**
	 * @brief Converts the letters, AABb, etc... into a numerical index
	 */
	uint GetGenotypeIdx(const char *genotype, uint lociCount);

	void ReportDiseaseLoci(ostream &os, vector<Locus *> &diseaseLoci);
	void ReportPenetranceTable(ostream &os);

	virtual void GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci) = 0;

	/**
	 * @brief Returns the number of loci associated with the model
	 */
	size_t GetModelSize() 							{ return loci.size(); }

	/**
	 * @brief Indicates that the model is continuous (or discrete)
	 */
	virtual bool IsContinuous() { return false;}

	/**
	 * @brief All penetrance based models can boil down to a simple look up for status
	 * @note more complex or non-boolean status will override this function
	 */
	int GetStatus(std::vector<uint>& genotypes, float& outcome);

	virtual bool Init(istream& s, PoolManager *pools)=0;

	virtual string GetModelConfiguration()=0;

	/**
	 * @brief initialize the penetrance at a given cell
	 */
	virtual void AddPenetrance(uint idx, double penetrance);
	DiseaseLocus GetLocus(istream& ss);
	/**
	 * @brief Helper factory method for parsing various configuration methods
	 * @param details Stream containing information about the model
	 * @param pools The pool manager which will provide additional details (al. Freq, etc)
	 * @return Penetrance based model (actual type is determined by details)
	 */
	static PenetranceModel *GetModel(istream& details, PoolManager *pools);
	

	static const char *PenTable;
	static const char *Simla;
	static const char *Simpen;
	static const char *Continuous;
	static const char *Template;

	void WritePenetranceFile(const char *filename, vector<Locus *>& diseaseLoci, float thresh);

	string ParseFilename(istream &s);

	int GetType() { return 1; }
	/**
	 * @brief Inverts the values associated with a given locus in a penetrance model

	static bool InvertAllele(vector<double>& table, int modelSize, int locusToInvert);
	 */	
protected:

	virtual int ConvertGenotype(int index, int gt);

	/**
	 * @brief Adds genotype labels into the vector (index of the label is the genotype value)
	 * @param labels The vector of labels associated with the set of genotypes
	 * @param modelSize The size of model associated with the penetrance table
	 */
	void BuildGenotypeLabels(vector<string> &labels, uint modelSize);

	/**
	 * @brief Converts an array of genotype values into the appropriate string
	 * @param genotypes Array of integers 0,1,2 representing the genotype at each locus
	 * @param number of values in the array above. 
	 */
	string BuildGenotypeLabel(uint *genotypes, uint modelSize);

	/**
	 * @brief Converts a single genotype to the appropriate letter based on it's position
	 * @param genotype 0,1,2 representing the genotype 
	 * @param which position (0 is first and will be encoded as A or a)
	 */
	string BuildGenotypeLabel(uint genotype, uint position);
	/**
	 * @brief sets up penetrance table details. Must be called once you know the number of loci
	 */
	void SetModelSize(size_t modelSize);

	/**
	 * @brief Used in determining the position with the model's penetrance table based on genotype
	 */
	uint GetMultiplier(uint genotype, uint position);

	/**
	 * @brief constructs the compressed genotype information used by simla
	 */
	int BuildGenotypeIndex(vector<uint> & genotypes);

	int BuildGenotypeIndex(int *genotypes, int numGenos);

	
	//For speed, we are using an array to represent these values
	uint penCount;						///<Size of the penList array
    double *penList;					///<The array of penetrances for each of the loci
	uint maxPenCount;					///<The maximum number of penetrances possible for a given sized model	


	bool isConfigured;					///<Determines if the model is properly configured
	uint observed;						///<Used to report observed prevalence
	uint totalQueries;					///<Total number of times IsAffected was called

};





inline
PenetranceModel::PenetranceModel() 
	: penCount(0), penList(NULL), maxPenCount(0), isConfigured(false), observed(0), totalQueries(0) { }




}


}
#endif 
