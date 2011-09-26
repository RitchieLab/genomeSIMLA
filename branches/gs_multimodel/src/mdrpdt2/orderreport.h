//
// C++ Interface: orderreport
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_EVALUATIONORDERREPORT_H
#define MDRPDT_EVALUATIONORDERREPORT_H

#include "utility/types.h"
#include "utility/rbtree.h"
#include "foldstatistic.h"
#include "ptestdistribution.h"

namespace MdrPDT {

namespace Evaluation {

class EvaluationMethod;

typedef Utility::RBTree<float, PdtModel, floatCompare> TReportTree;
typedef Utility::RBTreeNode<float, PdtModel, floatCompare> TReportTreeNode;

/**
@brief Manages reports for each order (useful for storing reports for each fold for each order...etc)

	@author Eric Torstenson
*/
class OrderReport {
public:
	OrderReport(int xvCount, int reportSize, int order);

	/**
	 * @brief Drop a new result into the report
	 */
	void AddModelResult(PdtModel& result);


	float GetAvgMOR(int idx);
	/**
	 * @brief Returns the best model at a given index
	 */
	FoldStatistic GetBestModel(int idx);

	/**
	 * @brief Generates a summary report for this order
	 */
	void GenerateSummaryReport(std::ostream& os, MdrPDT::Distribution::PTestDistribution* dist);

	/**
	 * @brief Generates a detailed for a given index into the report
	 */
	void GenerateDetailedReport(int idx, std::ostream& os, EvaluationMethod *eval);
	~OrderReport() { }

protected:
	std::vector<TReportTree> reportFolds;	///<This is the N records that represent the top runs
	int reportSize;							///<The number of records to keep
	int order;								///<Which order we are working with
};


}

}

#endif
