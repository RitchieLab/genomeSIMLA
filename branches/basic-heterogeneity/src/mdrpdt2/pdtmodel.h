/**
 * PdtModel 
 * 
 * Description:		This class encapsulates basic model details required for evaluating a model
 * 
 * Author: 		Eric Torstenson, (C) Marylyn Ritchie 2008
 * 
 * Copyright:		See COPYING file that comes with this distribution.
 */
#ifndef PDT_MODEL_H
#define PDT_MODEL_H

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "genotypecounter.h"
namespace MdrPDT {

using namespace std;




class PdtModel {
public:
	PdtModel();
	~PdtModel();
	
	/**
	 * @brief less that operator required for std sorting
	 */
	bool operator<(const PdtModel& other);
	
	/**
	 * @brief streams report details to os for a given fold
	 */
	void ReportFold(int id, std::ostream& os);
	
	/**
	 * @brief returns the training T score for the whole item
	 */
	float GetTrainingT(int id) const;
	/**
	 * @brief returns the testing T score for who item at id
	 */
	float GetTestingT(int id) const;
	
	/**
	 * @brief Evaluates a single fold, returning training & testing
	 */
	float EvaluateFold(float &training, float& testing, int id);
	
	/**
	 * @brief Adds a new snp to the model
	 */
	void AddSnp(int id, char *ptr);
	
	/**
	 * @brief Adds a new fold to the model
	 */
	void AddFold(Evaluation::tcalc& training, Evaluation::tcalc& testing);
	
	/**
	 * @brief returns a textual representation of the model
	 */
	std::string GetModelID();
	std::string GetPrintableLabel();
	/**
	 * @brief Returns the size of the model
	 */
	int ModelSize();
	
	/**
	 * @brief Returns the id for a given index (i.e. 5x10, the id for 1 is 10)
	 */
	int GetLocusID(int idx);
	
	/**
	 * @brief Returns the genotype information at genotype, idx
	 */
	char *GetSnpData(int idx);
	
	/**
	 * @brief Setup the pedigree details
	 */
	void InitPedigree(int indCount, int count, int *data);
	
	/**
	 * @brief Initialize the cross validation details
	 */
	void CrossValidations(int count, char *data);
	
	/**
	 * @brief Returns the array of folds associated with each individual
	 */
	char *FoldData();

	/**
	 * @brief Returns the array of pedigrees for each individual
	 */
	int *PedigreeData();
	
	/**
	 * @brief Returs the number of cross validation slices in use
	 */
	int FoldCount();
	
	/**
	 * @brief Returns the number of pedigrees associated with the data
	 */
	int PedigreeCount();
	
	/**
	 * @brief Number of individuals associated with the dataset
	 */
	int IndividualCount();
	
	/**
	 * @brief Returns the vector of Model IDs
	 */
	vector<int> GetModelIDs();
	
	/**
	 * @brief Returns the order of the model
	 */
	int GetModelSize();
	
	/**
	 * @brief Returns the ID assocaited with a given snp
	 */
	int GetSnpID(int n);

	string GetSnpLabel(int n);
	
	void ShowGenotypeIDs() {
		for (size_t i=0; i<loci.size(); i++) {
			cout<<"\t"<<i<<"\t\t";
			for (int o=0; o<individualCount; o++) {
				cout<<(int)loci[i][o]<<" ";
			}
			cout<<"\n";
		}
	}
	vector<char *> GetModelGenotypes() { return loci; }
	
	string GetGenotypeEncoding(int snpIdx, int genotype) {
		static string gts[] = {"1/1", "1/2", "2/2"};
		assert(genotype<=2 && genotype>=0);
		return gts[genotype];
	}
	
	void InitModelCounts();
	Evaluation::MatchedOddsRatio EvaluateMOR();
	Evaluation::MatchedOddsRatio EvaluateMOR(int fold);
	Evaluation::MatchedOddsRatio EvaluateBiasedMOR(int fold);
	Evaluation::GenotypeConverter& GetConverter() { return makeGenotype; }
	Evaluation::OverallCounter& GetGenotypeCounter() { return genotypeCounter; }
	void SetSnpLabels(std::vector<std::string>* labels) { snpLabels=labels; }
	void SetDSPCounts(std::vector<int>* counts) { dspCounts=counts; }
	int IndividualCount(int fold) { return (*dspCounts)[fold]; }

	static int MaxLabelLength;
protected:
	std::vector<char*>	loci;				///<The vector of genotype arrays
	string modelID;							///<model identification
	std::vector<Evaluation::tcalc> training;			///<training info for each fold
	std::vector<Evaluation::tcalc> testing;				///<Testing info for each fold
	std::vector<int>	snpIDs;				///<The indices associated with the model
	char *foldData;							///<array of fold IDs
	int *pedigreeData;						///<Array of pedigree IDs
	int xvCount;							///<The number of folds
	int pedCount;							///<The total number of pedigrees
	int individualCount;					///<The number of individuals that are part of the study
	Evaluation::GenotypeConverter makeGenotype;			///<Basic converter object for reporting (0->1/1 1->1/2 2->2/2)
	Evaluation::OverallCounter genotypeCounter;			///<Stash the counts for use by all evaluations
	
	std::vector<std::string> *snpLabels;
	std::vector<int> *dspCounts;
	
};

}

#endif
