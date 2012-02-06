//
// C++ Interface: configurationparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESECONFIGURATION_H
#define ESECONFIGURATION_H

#include "utility/utility.h"
#include "genetics/gtfileparserbuffered.h"
#include "gtlineparsermdrpdt.h"
#include "genetics/snpevalsuite.h"
#include "genetics/gtfileparser.h"
#include "reportlog.h"
#include "locuslogascii.h"
#include <string>

namespace MDRPDT {

using namespace std;
using namespace Utility;
using namespace Genetics::Parser;

/**
@brief Parses the text file which contains the various configuration details. 
Pass one of these objects to an instance of LineParser
@todo Figure out the probability of A and a (see configuration.cpp:82 from original ESE)

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class AppConfiguration : public ConfigurationParser {
public:
	/**
	 * @brief Used to record the style of analysis is to be performed
	 */
	typedef enum AnalysisStyle{
		NoAnalysis, 						///<This will eventually used if you need to perfrom some basic processing or load a pre-analyized dataset
		MaxDiff,							///<Perform Basic ESE evaluation- max diff
		BalancedAccuracy,
		PDT
	} AnalysisStyle;
	
	/**
	 * @brief Used to record the format the input file is in
	 */
	typedef enum {
		InputSpaceDelimited,
		InputLabeledSpaceDelimited,
		InputBinaryInput,
		MdrFormat,
		Pedigree
	} InputFormat;


    AppConfiguration();						///<Constructor
    ~AppConfiguration();					///<Destructor
	/**
	 * @brief return the file parser associated with the current setup
	 */
	GtFileParser *GetInputFileParser();

	/**
	 * @brief Returns the PDT specialized evaluation suite. 
	 * @param doRandomizeStatus Set for true if the evaluation is for a permutation test
	 * @param testNumber Used to indicate which test we are currently on
	 */
	SnpEvalSuite    *GetPdtEval(bool doRandomizeStatus, uint testNumber);

	/**
	 * @brief Allow the application to close the parser so that we can reclaim any buffered memory, if needed
	 */
	void CloseInputParser();			
	
	/**
	 * @brief Indicates that we are to perform a perfect search (or not)
	 */
	bool DoPerfectSearch() { return performPerfectSearch; }

	/**
 	 * Dumps a Text Report summarizing the configuration options encountered
	 */
	void GenerateReport(ostream& os);	

	/**
	 * Setup the key value assignments
	 */
	bool SetValues(const char *key, const char *value);

	/**
	 * @brief Perform tests to make sure everything that is set in the configuration makes sense
	 */
	bool Validate();


	/**
	 * @brief Performs any post loading work that needs to be done
	 * Always call this after configuration file has been loaded
	 */
	void PostLoad();

	static AppConfiguration *Instance();		///<Singleton reference acquisition
	static void Release();						///<Singleton Release

	/**
	 * @brief Return the Distribution log
	 */
	BasicLog *GetDistributionLog();			
	BasicLog *GetReportLog();					///<Return the report log (used for writing all analysis results)
	BasicLog *GetPedigreeLog();					///<Used to log any pedigree related information

	//ESE Requirements
	uint affectedCount	;						///<The number of affected individuals
	uint unaffectedCount;						///<The number of unaffected individuals
	uint bestDifferenceThreshold;				///<The threshold for selecting models of interest
	uint reportModelCount;						///<The number of models to be added to the final report (0 will report all)

	//General settings
	string inputFile;							///<Where the snps come from	- Later, this will be a list
	bool doSortReport;							///<Indicate that the model report should be sorted (only checked if all models are reported)
	uint pTestCount;							///<The number of PTests to be performed
	float pValueThreshold;						///<Threshold used for reporting with p-values
	float minDThresh;							///<Threshold for the D range optimization
	uint pTestSeed;								///<The seed used for PTests

	AnalysisStyle analysisStyle;				///<Style of analysis required
	LocusLogAscii::LocusWriteStyle locusWriteStyle;
	
	LocusLogAscii *logLocus;					///<The locus log (if being used)
	InputFormat inputFormat;					///<The format the input file is in

	string errorMsg;							///<Used to report problems encountered in certain translations from config file

	uint comboEnd;								///<The largest size of models to be considered
	uint comboStart;							///<The minimum size of models to be considered

	uint maxGenotypes;							///<The maximum number of gennotypes to be found
	PermutationTestDist *distribution;			///<This is the distribution that is being used


	/**
	 * @brief Singleton reference counting
	 */
	static uint _instanceCount;

	/**
	 * @brief Singleton pointer
	 */
	static AppConfiguration *_instance;
	

	/**
	 * @brief The family repository
	 */
	FamilyRepository *families;

	/**	
 	 * @brief The current evaluation suite
	 */
	SnpEvalSuite *evalSuite;
	bool doCreateVirtualSibs;

	string reportName;							///<Prefix used for the all reports
	bool doCalcThreshold;
	struct ConfigFileValues {
	public:
		//Fundamental settings
		static char *ComboStart;				///<The first locus count for models of interest
		static char *ComboEnd;					///<The largest size of interest
		static char *CalcThreshold;				///<Instructs the application to evaluate risk based on model threshold

		//Settings for overrideing various report extensions
		static char *ExtDistribution;			///<Let user set the extension for distribution reports
		static char *ExtReport;					///<Let user set the extension for the overall report
	 	static char *ExtLocus;					///<Let user set the extension used for locus reports
		static char *ExtPedigree;				///<Let user set the extension used for the pedigree report
		static char *ReportName;				///<Let user set the report name
		
		//Overriding other reporting details
		static char *MinPVal;					///<Let the user assign the P-Value cutoff for reporting
		static char *PTestCount;				///<Number of PTests
		static char *PTestSeed;					///<Seed associated with the tests
	
		//Overriding input details
		static char *_InputFile;				///<Used to override the default input format
		static char *_InputFormat;				///<Used to specify the format needed for the input file
		static char *_MdrFormat;				///<Used to specify to use MDR formatted input
		static char *PedigreeFormat;			///<Used to specify the pedigree format

		//Overriding how the data is to be interpretted
		static char *_GenotypeCount;			///<Use to override the max genotypes
		static char *_AnalysisStyle;			///<Used to override the type of analysis
		static char *BalancedAccuracy;			///<Used to specify Balanced Accuracy Analyses on CC data
		static char *PDTAnalysis;				///<Used to specify PDT analyses on pedigree data
		static char *_MaxDiff;					///<Used to specify Max Difference on CC data
		
		//MD Specific Settings
		static char *_ReportModelCount;			///<Used to overide the model count desired for final report
		static char *_AffectedCount;			///<Used to specify the number of affected individuals (in the Base Pair format)
		static char *_UnaffectedCount;			///<Used to specify the nubmer of unaffected individuals (in the base pair format)
		static char *_BestDifferenceThreshold;	///<Threshold used to determine when a model is of interest
		static char *_SortReportModels;			///<Indicates to sort the reported models
		
		//Pedigree Specific settings
		static char *MinD;						///<The percentage within the highest D to be considered
		static char *CreateVirtualSibs;			///<Toggles whether virtual sibs are created for triads

		static char *AffectedValue;				///<Used to specify the affected value
		static char *UnaffectedValue;			///<Used to specify the unaffected value

		static char *VerboseMOR;				///<Verbose Matched Odds Ratio
		static char *VerbosePDT;				///<Produces the statistic for each model evaluated

	};		

protected:

	void SetAnalysisStyle(const char *value);	///<Used to set the analysis style indicated by the configuration file
	void SetInputFormat(const char *value);		///<Used to set input format indicated by the configuration file
	void SetDefaultValues();					///<Used to set up the default values
	GtFileParser *inputParser;
	string extDistribution;						///<Extension used to set up distribution report filename
	string extReport;							///<Extension used for the overall report
	string extLocus;							///<Extension used for locus reports
	string extPedigree;							///<Extension used for the pedigree report

	//Wrappers for some of the logs
	BasicLog *logDistribution;
	BasicLog *logReport;
	BasicLog *logPedigree;
	bool performPerfectSearch;					///Used to specify if perfect search is desired

};
inline
void AppConfiguration::CloseInputParser() {
	if (inputParser)
		delete inputParser;
	inputParser=NULL;
}

inline
AppConfiguration *AppConfiguration::Instance() {
	if (_instance == NULL)
		_instance = new AppConfiguration();
	_instanceCount++;
	return _instance;
}

inline
void AppConfiguration::Release() {
	if (--_instanceCount < 1) {
		delete _instance;
		_instance = NULL;
	}
}

}

#endif