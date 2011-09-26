//
// C++ Implementation: pedigreestatistics
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pedigreestatistics.h"
#include <iomanip>

using namespace std;

namespace MdrPDT {

namespace Validator {

PedigreeStatistics::PedigreeStatistics(ostream& summary)	: 
		summaryReport(summary) { }


PedigreeStatistics::~PedigreeStatistics()	{
}


/**
 * @brief Perform tasks associated with a specific individual
 */
void PedigreeStatistics::Evaluate(Pedigree* ped, Individual *ind) {
	if (ped->DropFromAnalysis() || ind->DropFromAnalysis()) {
		PedigreeContributions& c=nonParticipating[ped->ID()];
		if (ind->IsAffected()) {
			c.individualCounts.first++;
			totalIndividuals.individualCounts.first++;
		} else {
			c.individualCounts.second++;
			totalIndividuals.individualCounts.second++;
		}
	}
	else {
		PedigreeContributions& c=participating[ped->ID()];
		if (ind->GetMother() == NULL || ind->GetFather() == NULL) {
			if (ind->IsAffected()) {
				c.individualCounts.first++;
				founders.individualCounts.first++;
				totalIndividuals.individualCounts.first++;
			}
			else {
				c.individualCounts.second++;
				founders.individualCounts.second++;
				totalIndividuals.individualCounts.second++;
			}
		}		
		else if (ind->IsAffected()) {
			c.individualCounts.first++;
			participatingIndividuals.individualCounts.first++;
			totalIndividuals.individualCounts.first++;
		} else {
			c.individualCounts.second++;
			participatingIndividuals.individualCounts.second++;
			totalIndividuals.individualCounts.second++;
		}
	}
}
	
/**
 * @brief Prepare to start an evaluation (basic book keeping)
 */
void PedigreeStatistics::PrepEvaluation() {
	participating.clear();
	nonParticipating.clear();
	participatingIndividuals.Reset();
	totalIndividuals.Reset();
}

/**
 * @brief Finish up with the details of an evaluation
 */
void PedigreeStatistics::PostEvaluation() {
	//Summarize the counts
	summaryReport<<setw(45)<<"Total Individuals: "
			<<totalIndividuals.individualCounts.first+totalIndividuals.individualCounts.second<<"\n";
	summaryReport<<setw(45)<<"Participating Individuals: "
			<<participatingIndividuals.individualCounts.first+participatingIndividuals.individualCounts.second<<" ( "
			<<participatingIndividuals.individualCounts.first<<"A | "
			<<participatingIndividuals.individualCounts.second<<"U )\n";
	summaryReport<<setw(45)<<"Participating Founders: "
			<<founders.individualCounts.first + founders.individualCounts.second<<"\n";


}

}

}
