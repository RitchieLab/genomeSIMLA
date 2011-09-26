//
// C++ Interface: powerstudy
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPOWERSTUDY_H
#define MDRPOWERSTUDY_H
#include "snpsearchapplication.h"
#include "eseconfiguration.h"
#include "genetics/snprepository.h"
#include <fstream>

namespace MDR {

namespace Power {
/**
@brief Application designed to perform power studies on MDR-PDT

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class PowerStudy : public SnpSearchApplication {
public:

struct ModelResults {
	ModelResults() : freq(0), sumCls(0.0), sumPrd(0.0), sumTraining(0.0), sumTesting(0.0), sumOR(0.0), totalT(0.0),
	trainingFreq(0), testingFreq(0), clsFreq(0), prdFreq(0), orFreq(0) { }

	double GetAverageT() {
		return totalT / (float)freq;
	}

	ModelResults &operator+=(ModelStatistics& other) {
		freq++;
		if (!isnan(other.GetAvgClassificationError())) {
			sumCls+=other.GetAvgClassificationError();
			clsFreq++;
		}
		if (!isnan(other.GetAvgPredictionError())) {
			sumPrd+=other.GetAvgPredictionError();
			prdFreq++;
		}
		if (!isnan(other.GetAvgTesting())) {
			sumTesting+=other.GetAvgTesting();	
			testingFreq++;
		}
		if (!isnan(other.GetAvgTraining())) {
			sumTraining+=other.GetAvgTraining();
			trainingFreq++;
		}
		if (!isnan(other.GetOddsRatio())) {
			sumOR+=other.GetOddsRatio();
			orFreq++;
		}
		return *this;
	}

	float GetAverageTraining() { 	return sumTraining / (float)trainingFreq;	}
	float GetAverageTesting()  {	return sumTesting / (float)testingFreq; }
	float GetAverageCls() { 		return sumCls / (float)clsFreq; }
	float GetAveragePrd() { 		return sumPrd / (float)prdFreq; }
	float GetAverageOR()  {			 return sumOR  / (float)orFreq; }

	uint freq;
	uint foldCount;
	float sumCls;
	float sumPrd;
	float sumTraining;
	float sumTesting;
	float sumOR;
	float totalT;

	uint trainingFreq, testingFreq, clsFreq, prdFreq, orFreq;
	
};

struct ModelHolder {
	ModelHolder(const char *filename, uint pos, uint size, ModelStatistics &st) : file(filename), stat(st), position(pos), modelSize(size) {}
	~ModelHolder() {}

	void GenerateHeader(ostream* os) {
		*os<<setw(29)<<right<<"File Name"<<" ";
		*os<<setw(8)<<left<<"Order";
		*os<<setw(15)<<"Model ID";
		if (stat.foldCount > 1) {
			*os<<right<<setw(10)<<"XV Cons.";
			*os<<right<<setw(10)<<"Cl. Error";
			*os<<right<<setw(10)<<"Training";
			*os<<right<<setw(10)<<"Pr. Error";
			if (pDist)
				*os<<right<<setw(10)<<"Pr. P-Value";
			*os<<right<<setw(10)<<"Testing";
			if (fDist)
				*os<<right<<setw(10)<<"T. P-Value";
		}
		else {
			*os<<right<<setw(10)<<"Cl. Error";
			*os<<right<<setw(10)<<"MDR-PDT";
			if (fDist)
				*os<<right<<setw(10)<<"T. P-Value";
		}
		*os<<right<<setw(10)<<"MOR";
		if (oDist)
			*os<<right<<setw(10)<<"MOR PV";
		*os<<"\n";
	}
	void GenerateReport(ostream* os) {
		*os<<setw(29)<<right<<file<<" ";
		*os<<setw(8)<<left<<modelSize + 1;
		*os<<setw(15)<<stat.label;
		if (stat.foldCount > 1) {
			*os<<right<<setw(10)<<stat.xvConsistency;
			*os<<right<<setw(9)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<stat.GetAvgClassificationError()<<"%";
			*os<<right<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<stat.GetAvgTraining();
			*os<<right<<setw(9)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<stat.GetAvgPredictionError()<<"%";
			if (pDist)
				*os<<right<<setw(10)<<setprecision(4)<<pDist->GetPValue(stat.GetAvgPredictionError(), modelSize);

			*os<<right<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<stat.GetAvgTesting();
			if (fDist)
				*os<<right<<setw(10)<<setprecision(4)<<fDist->GetPValue(stat.GetAvgTesting(), modelSize);

		}
		else {
			*os<<right<<setw(9)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(2)<<stat.GetAvgPredictionError()<<"%";
			*os<<right<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<stat.GetAvgTesting();
			if (fDist)
				*os<<right<<setw(10)<<setprecision(4)<<fDist->GetPValue(stat.GetAvgTesting(), modelSize);
		}
		*os<<right<<setw(10)<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<stat.GetOddsRatio();
		if (oDist)
			*os<<right<<setw(10)<<setprecision(4)<<oDist->GetPValue(stat.GetOddsRatio(), modelSize);

		*os<<"\n";
	}

			
	
	string file;						///<This is the filename
	ModelStatistics stat;				///<The leading statistic for this one
	uint position;						///<If we are reporting more than the top models, this is used to differentiate
	uint modelSize;						///<The order of models where this occured

	/**
	 * Copies of the distributions. 
	 */
	static PTestDistribution 		*fDist;
	static PTestDistribution 		*oDist;
	static PTestDistribution 		*pDist;

};


    PowerStudy();

    ~PowerStudy();



	/**
 	 * @brief This will let the program set up the necessary parameters based on cmd line arguments. 
	 * @return true indicates that execution can begin
 	 */
	virtual bool ParseCmdLine(int argc, char** argv);

	/**
	 * @brief Prints the help contents
	 */
	virtual void PrintHelp();
	
	/**	
	 * @brief Starts execution
	 */
	virtual void Start();
	
	/**
	 * @brief Threadsafe way to determine which dataset is coming next
	 * @note Serial & MPI/Master nodes simply read next file from the array
	 * @note MPI/Salve nodes wait for message from Master node to determine which to use
	 */ 
	bool NextDataset(string& filename);

	/**
	 * @brief Performs various tasks associated with ref. distribution based on node type
	 * @note For MPI/Master nodes, this will save the distribution to files, wait for slaves to finish loading
	 * 		and delete the files
	 * @note For MPI/Slave nodes, the nodes will wait for an array of filenames to load the various distributions
	 * @note For Serial, nothing happens
	 */
	void ManageReferenceDist();

	static bool VerbosePower;

#ifdef USE_MPI
	/**
	 * @brief Initiates the main thread to manage the various assignments to slave nodes (then sleep). 
	 */
	static void *ProcessWatcher(void *arg);

	/**
	 * @brief Initiates thread to set up the reference distribution.
	 */
	static void *PTestWatcher(void *arg);
#endif

protected:
	/**
	 * @brief Saves the result of a run in the appropriate way (serial, MPI/Master, MPI/Savel)
	 */
	void SaveResults(uint nodeID, const char *filename, uint position, uint modelSize, ModelStatistics &st);

	/**
	 * @brief Responsible for sending a p-test to the Master Node
	 */
	void SavePTest(int nodeID, int testID, int modelSize, ModelStatistics &st);

	/**
	 * @brief Responsible for aggregating a complete run into the various distributions
	 * @note Fitness Distribution is N-Test, so each order model is treated to it's own distribution
	 * @note MOR / PE distributions are based on omnibus. The winning model is based on the following, in order:
	 * 			Cross Validation Consistency, MOR / PE value (whichever distribution we are building)
	 * @note For Serial & MPI/Master, simply add the new node to the distribution
	 * @note For MPI/Slave, the values are sent to the Master
	 */
	void SavePTests(int nodeID, int testID, vector<ModelStatistics> bestModels);
	/**
	 * @brief Do a quick sanity check on the settings provided
	 */
	bool VerifyConfiguration();
	string targetModel;
	vector<string> inputfiles;

	vector<ModelHolder> statistics;				///<The statistics for each run. These are in order
	map<string, ModelResults> power;			///<Overall power of the fitness evaluation
	map<string, ModelResults> orPower;			///<Power using Matched Odds Ratio
	map<string, ModelResults> pePower;			///<Power based on prediction error
	uint fileIdx;								///<Used to keep up with the index of the current file
	string configFile;

	bool useReferenceDistribution;				///<When true, we will build one set of distributions and each dataset will
												///<reference that distribution 

	/**
	 * @brief Basic pedigree search over a single dataset
	 */
	double PedigreeSearch(const char *filename, uint seed, uint testID);

	/**
	 * @brief Handles the reporting of a single power result structure
	 * @note Dumps the contents of the map as well as the power/details of the target model
	 */
	void ReportPower(map<string, ModelResults>& power);


};

}
}
#endif
