//
// C++ Interface: locusreport
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLOCUSREPORT_H
#define SIMULATION_VISUALIZATIONLOCUSREPORT_H
#include "locus.h"
#include "blocklistnode.h"
#include "locselection.h"

namespace Simulation {

namespace Visualization {


/**
@ brief Providing a simple and clear representation of the loci of interest

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LocusReport{
public:
    LocusReport();	
	LocusReport(const LocusReport& other) : selectors(other.selectors) { }
	
    ~LocusReport();
 
	void AddLocusSelection(LocusSelection &lc);

	/**
	 * @brief Initialize the report with the loci as well as identify related blocks they are present in
	 */
	void Init(vector<Locus>& loci, BlockList *blocksOfConcern, const char *reportFilename, const char *chromID);

	/**
	 * @brief Writes html report to the stream
	 * @param file This is the stream to be written to
	 * @param sisterFile The name of the file containing the block reports (links will be made to this file)
	 */
	void WriteHtmlReport(ostream &file);

	
	void Reset();
protected:
	vector<LocusSelection> selectors;
};


inline
LocusReport::LocusReport() { }

inline
void LocusReport::AddLocusSelection(LocusSelection &lc) {
	selectors.push_back(lc);
}




}

}

#endif
