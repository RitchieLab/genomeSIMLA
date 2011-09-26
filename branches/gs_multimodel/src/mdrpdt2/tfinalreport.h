//
// C++ Interface: tfinalreport
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_EVALUATIONTFINALREPORT_H
#define MDRPDT_EVALUATIONTFINALREPORT_H
#include <string>
#include "evaluationmethod.h"
#include "orderreport.h"
#include "ptestdistribution.h"

class EvalutionMethod;
namespace MdrPDT {

namespace Evaluation {

/**
@Brief Report for T-Based results

	@author Eric Torstenson
	@brief Aggregates the various PdtModels and keeps the best N for each fold
*/
class TFinalReport {
public:
	TFinalReport(int maxModelSize, int xvCount, int reportSize) ;
	~TFinalReport();

	/**
	 * @brief Adds a the result of an evaluation to the report		
	 * @Note By default, the on ReportSize results will be saved
	 */
	void AddModelResult(PdtModel& result);

	/**
	 * @brief Return the top MOR for a given order
	 */
	float GetBestMOR(int order, std::string& id);
	float GetBestMOR(std::string& id);
	float SetAvgMOR(float mor);
	/**
	 * @brief generate verbose reports. 
	 * @param os the stream to be written to
	 * @param depthOfReport The number of reports to produce (1 is top only)
	 */
	void GenerateVerboseReport(std::ostream& os, int depthOfReport, EvaluationMethod *eval);

	/**
	 * @brief Generates basic resport based on each entry in the report
	 */
	void GenerateSummaryReport(std::ostream& os, Distribution::PTestDistribution* dist);

	/**
	 * @brief Generates verbose report for the model at index idx where 0 has the highest testing score
	 */
	void GenerateDetailedReport(int idx, std::ostream& os, EvaluationMethod* eval);

	static bool VerboseAnalysis;
protected:
	std::vector<OrderReport> reports;			///<The order reports for each order
	int reportSize;								///<The number of reports kept
};


}

}

#endif
