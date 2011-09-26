//
// C++ Interface: pedigreestatistics
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MDRPDT_VALIDATORPEDIGREESTATISTICS_H
#define MDRPDT_VALIDATORPEDIGREESTATISTICS_H
#include <iostream>
#include <map>
#include "pedigreerepositorytraverser.h"

namespace MdrPDT {

namespace Validator {

/**
@brief Gathers basic statistics and reports on a repository's participants at both individual and pedigree level

	@author Eric torstenson
*/
class PedigreeStatistics : public PedigreeRepositoryTraverser {
public:
    PedigreeStatistics(std::ostream& summary);
    ~PedigreeStatistics();

	/**
	 * @brief Perform tasks associated with a specific pedigree
	 * @note There is no need to call Evaluate(Individual), that is done by another function
	 */
	void Evaluate(Pedigree* ped, Individual* ind);
	
	/**
	 * @brief Prepare to start an evaluation (basic book keeping)
	 */
	void PrepEvaluation();

	/**
	 * @brief Finish up with the details of an evaluation
	 */
	void PostEvaluation();
	


	/**
	 * @brief 
	 * @note First in each pair is affected
	 */
	struct PedigreeContributions {
		std::pair<int,int> individualCounts;
		std::string id;

		PedigreeContributions(const char *id="", int aff=0, int unaff=0) : 
			individualCounts(aff, unaff), id(id){ }

		void Reset() {
			individualCounts=std::pair<int,int>(0,0);
			id="";
		}
		PedigreeContributions& operator+(const PedigreeContributions& other) {
			individualCounts.first+=other.individualCounts.first;
			individualCounts.second+=other.individualCounts.second;
			id += " " + other.id;	
			return *this;
		}
	};
protected:
	std::ostream& summaryReport;
	PedigreeContributions founders;
	///Keep up with each participating pedigree members
	std::map<std::string, PedigreeContributions> participating;
	///Keep up with those not participating
	std::map<std::string, PedigreeContributions> nonParticipating;
	///Just sum over all...for convenience
	PedigreeContributions totalIndividuals;
	///Just sum over all contributing pedigrees
	PedigreeContributions participatingIndividuals;
};


}

}

#endif
