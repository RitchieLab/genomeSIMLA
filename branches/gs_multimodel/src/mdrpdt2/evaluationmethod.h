//
// C++ Interface: evaluationmethod
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_EVALUATIONEVALUATIONMETHOD_H
#define MDRPDT_EVALUATIONEVALUATIONMETHOD_H
#include <iostream>
#include "genotyperepository.h"
#include "pdtmodel.h"
namespace MdrPDT {

namespace Evaluation {

class TFinalReport;
/**
@brief Baseclass for evaluation methods

	@author Eric Torstenson
*/

class EvaluationMethod {
public:
	EvaluationMethod() { }
	virtual ~EvaluationMethod() { }


	virtual void EvaluateModel(PdtModel& report)=0;
	virtual void EvaluateModelVerbose(PdtModel& report, ostream& os) = 0;
	virtual void EvaluateModelVerbose(PdtModel& report, int fold, ostream& os) = 0;
	virtual void EvaluateModel(GenotypeRepository *rep, const char *model)=0;
	/**
	 * @brief Perform evaluation over all of the SNPs in the repository
	 */
	virtual void BasicEval(GenotypeRepository *rep, TFinalReport& report) = 0;
	virtual void BasicEval(GenotypeRepository *rep, TFinalReport& report, std::vector<char *>snps, std::vector<int> modelIDs, int snpID) = 0;

//	virtual void BasicEval(int *pedID, char *folds, char *snp1, PdtReport& report) = 0;


	/**
	 * @brief Perform verbose evaluation for reporting purposes
	 */
//	virtual void VerboseEval(int *pedID, char *folds, char *snp, std::ostream& os) = 0;

	static int maxModelSize;
	static int minModelSize;
};

}

}

#endif
