//
// C++ Interface: tstastic
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDTTSTASTIC_H
#define MDRPDTTSTASTIC_H

#include <deque>
#include "utility/types.h"
#include "evaluationmethod.h"
#include "tfinalreport.h"
namespace MdrPDT {

namespace Evaluation {

/**
@brief Evaluation method(s) for MDR-PDT

	@author 
*/
class TStatistic : public EvaluationMethod {
public:
    TStatistic(int xvCount, int pedCount, int indCount);
    ~TStatistic();
	
	void EvaluateModel(PdtModel& model);

	void EvaluateModelVerbose(PdtModel& report, ostream& os);

	void EvaluateModelVerbose(PdtModel& report, int fold, ostream& os);
	/**
	 * @brief Iterate over each of the various combinations and record their results
	 */
	void BasicEval(GenotypeRepository *rep, TFinalReport& report );


	/**
	 * @brief Constructs a single model based on the string, model, and performs a verbose evaluation
	 */
	void EvaluateModel(GenotypeRepository *rep, const char *model);
	/**
	 * @brief Evaluate a single model and return it's value
	 */
//	void BasicEval(int *pedID, char *folds, char *snp1, PdtReport& report);

	void BasicEval(GenotypeRepository *rep, TFinalReport& report, std::vector<char *>snps, std::vector<int> modelIDs, int snpID);

	/**
	 * @brief Evaluate a single model, but report each step to the stream for reporting purposes

	void VerboseEval(int *pedID, char *folds, char *snp, std::ostream& os);
 	 */	
protected:
	void GetGenotype(int genotype, int locusCount, std::deque<int>& genotypes);
	int xvCount;
	int pedCount;
	int indCount;


	int snpCount;
	char *folds;
	int *pedIDs;
};

}

}

#endif
