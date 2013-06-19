//
// C++ Interface: diseasemodeldetails
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENOMESIM_GUIDISEASEMODELDETAILS_H
#define GENOMESIM_GUIDISEASEMODELDETAILS_H

#include "utility/lineparser.h"
#include "simulation/locus.h"
namespace GenomeSIM {

namespace GUI {

using namespace Simulation;
using namespace Utility;

/**
	@brief Base class for importing details about disease models from files
	@note This is a more mature approach to importing data from the file (assuming the
			file isn't terribly large.) I will eventually replace the others with these 
			classes as time permits. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PenetranceModelDetails;
class SimlaModelDetails;
class SimpenModelDetails;

class DiseaseModelDetails{
public:
    DiseaseModelDetails(const char *filename);
    virtual ~DiseaseModelDetails();

	/**	
 	 * @brief Just populates the map with the details of the file and calls PostLoad
	 */
	virtual void Load();


	/**
	 * @brief write out each of the keywords and their data to the file
	 */
	virtual void Save(const char *filename) = 0;
	virtual void Save();

	/**
	 * @brief Generates the configuration details necessary for compilation. 
	 */
	virtual string GetConfigurationDetails(vector<Locus*>& loci)=0;
	
	/**
	 * @brief Open a model based on it's extension
	 */
	static DiseaseModelDetails *OpenDiseaseModel(const char *filename);

	/**
	 * @brief Open a model based on internal knowledge of what kind it should be
	 */
	static DiseaseModelDetails *OpenDiseaseModel(const char *type, const char *filename);

	/**
	 * @brief Indicate the allele frequency range for allele 1 on given locus
	 */
	virtual string GetLocusRange(int locIdx);
	virtual bool GetLocusRange(int locIdx, float &min, float &max, float &goal)=0;


	/**
	 * @brief Performs a basic sanity check (like ridiculous penetrances). 
	 * @param error (only set when return is false) contains details about the error
	 * @note This is not based on locus details
	 */
	virtual bool ValidateCfg(string& error) { return true; }
	/**
	 * @Brief The type of model
	 * @return PENTABLE, SIMPEN, SIMLA
 	 */
	virtual string GetType() = 0;

	/**
	 * @brief Returns the number of loci associated with a given disease. 
	 * @return -1 if the locus count is variable. Otherwise, it's the value returned
	 */
	virtual int GetLocusCount() = 0;
protected:
	/**
	 * @brief data will be fully populated and each subclass will extract relevant details
	 */
	virtual void PostLoad(FileToMap& data)=0;

	string filename;
};

class PenetranceModelDetails : public DiseaseModelDetails {
public:
	struct Freqs {
		float al1;
		float al2;
		char label;
		Freqs() : al1(0.5), al2(0.5), label('A') { }
		Freqs(float al1, float al2, char label) : al1(al1), al2(al2), label(label) { }
	};
	struct Penetrance {
		string label;
		double pen;
		Penetrance() : label(""), pen(0.0) { }
		Penetrance(const char *label) : label(label), pen(0.0) { }
		Penetrance(const char *label, double pen) : label(label), pen(pen) { }
	};
	PenetranceModelDetails(const char *filename);
	~PenetranceModelDetails();

	float freqThreshold;
	string GetConfigurationDetails(vector<Locus*>& loci);

	string GetType() { return "PENTABLE";}
	int GetLocusCount() { return locusCount; }

	bool ValidateCfg(string& error) {
		float sum = 0.0;
		int penCount = penetrances.size();
		if (penCount == 0) {
			error = "There are no cells associated with this penetrance. Please make sure the file is correctly configured.\n";
			return false;
		}
		for (int i=0; i<penCount; i++) {
			sum+=penetrances[i].pen;
		}
		if (sum == 0.0 || sum == penCount) {
			error = "The selected penetrance table appears to be designed in such a way as to make it impossible to produce either affecteds or unaffecteds. Please examine this file and correct the cells.";
			return false;
		}
		return true;
	}
protected:
	string BuildGenotypeLabel(uint genotype, uint position);
	string BuildGenotypeLabel(uint *genotypes);
	void BuildGenotypeLabels();
	void PostLoad(FileToMap& data);
	void LoadLoci(vector<string>& freqs);
	void LoadPenetrances(FileToMap& data);
	void Save(const char *filename);
	bool GetLocusRange(int locIdx, float &min, float &max, float &goal);
	vector<Freqs> frequencies;
	int locusCount;	
	vector<Penetrance> penetrances;
};




class SimpenModelDetails:public DiseaseModelDetails {
public:
	SimpenModelDetails(const char *filename);
	~SimpenModelDetails();
	void PostLoad(FileToMap& data);
	void Save(const char *filename);


	double heritWeight;
	double heritTarget;
	
	double varianceWeight;
	double varianceTarget;
	
	double orWeight;
	double orTarget;

	double penTarget;
	string gaSettings;
	string GetType() { return "SIMPEN";}
	int GetLocusCount() { return -1; }
protected:
	void WriteToFile(ostream& file, const char *key, vector<string>& values);
	string WriteSimpenConfigFile();
	string GetConfigurationDetails(vector<Locus*>& modelLoci);
	/**
	 * @brief this doesn't matter for simpen models. 
	 */
	bool GetLocusRange(int locIdx, float &min, float &max, float &goal) { return false; }
};


class SimlaModelDetails: public DiseaseModelDetails {
public:
	struct LocusOR {
		char label;
		bool diseaseAtMajor;
		double oddsRatio;
		double modelType;
		LocusOR() : label('A'), diseaseAtMajor(false), oddsRatio(0.0), modelType(0.0) { }
		LocusOR(char label, bool diseaseAtMajor, double oddsRatio, double modelType) : 
			label(label), diseaseAtMajor(diseaseAtMajor), oddsRatio(oddsRatio), modelType(modelType) { }
	};
	struct Interaction {
		string label;
		double oddsRatio;
		Interaction() : label(""), oddsRatio(0.0) { }
		Interaction(const char *label, double &oddsRatio) : label(label), oddsRatio(oddsRatio)  {} 
		string ConvertToSIMLA() {
			char letter = 'A';
			char num = '1';

			string newLabel(label);
			for (int i=0; i<10; i++) {
				size_t pos = newLabel.find(letter, 0);
				if (pos != string::npos) {
					newLabel = newLabel.replace(pos, 1, 1, num);
				}
				letter++;
				num++;
			}
			return newLabel;
		}
	};

	SimlaModelDetails(const char *filename);
	~SimlaModelDetails();
	void PostLoad(FileToMap& data);
	void Save(const char *filename);

	int locusCount;
	float prevalence;

	vector<LocusOR> 	loci;
	vector<Interaction> interactions;

	double ToBeta(double value) { 		return log(value);	}
	double FromBeta(double value) { 	return exp(value); }
	string GetType() { return "SIMLA";}
	int GetLocusCount() {  return loci.size(); }
protected:
	string GetConfigurationDetails(vector<Locus*>& modelLoci);
	bool GetLocusRange(int locIdx, float &min, float &max, float &goal) { return false; }

};


inline
string DiseaseModelDetails::GetLocusRange(int locIdx) {
	float min, max, goal;
	stringstream ss;
	if (GetLocusRange(locIdx, min, max, goal))
		ss<<goal<<" ("<<min<<" - "<<max<<")";
	else
		ss<<"Any";
	return ss.str();
}		
inline
bool PenetranceModelDetails::GetLocusRange(int locIdx, float &min, float &max, float &goal) {
	bool success=locIdx < frequencies.size();
	if (success) {
		goal=frequencies[locIdx].al1;	
		min=goal-freqThreshold;
		if (min < 0.0)
			min=0.0;
		max=goal+freqThreshold;
		if (max>1.0)
			max=1.0;
	}
	return success;
}		
	

}

}

#endif
