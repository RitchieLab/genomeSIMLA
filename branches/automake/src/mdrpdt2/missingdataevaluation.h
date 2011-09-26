//
// C++ Interface: missingdataevaluation
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_EVALUATIONMISSINGDATAEVALUATION_H
#define MDRPDT_EVALUATIONMISSINGDATAEVALUATION_H

#include <iostream>
#include "genotyperepository.h"

namespace MdrPDT {

namespace Evaluation {

/**
@Brief Generates a report on missing data for the repository. Also drops SNPs from analysis if they exceed the threshold for valid markers

	@author 
*/
class MissingDataEvaluation{
public:
    MissingDataEvaluation(float threshold);
    ~MissingDataEvaluation();

	/**
	 * @brief Perform analysis and generate report on the contents of the repository
	 */
	int Analyze(std::ostream& os, GenotypeRepository& repo);

protected:
	float threshold;
};

}

}

#endif
