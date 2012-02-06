//
// C++ Interface: gtfileparserbuffered
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_PARSERGTFILEPARSERBUFFERED_H
#define GENETICS_PARSERGTFILEPARSERBUFFERED_H

#include "gtfileparser.h"
#include "gtlineparser.h"
#include <iomanip>

namespace Genetics {

namespace Parser {



/**
@brief Opens up the data file(s) and fully loads all snp data during open(). Calls to NextSnp() returns a snp from memory. 

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GtFileParserBuffered : public GtFileParser {
public:
    GtFileParserBuffered(GtLineParser *parser, const char *filename);

    ~GtFileParserBuffered();

	/**
	 * @brief Returns the next snp
	 * @return SnpAligned pointer
	 */
	SnpAligned *NextSnp();
	/**
	 * @brief Opens up each file in the current group
 	 * @return True/False depending on wither there were problems opening the file(s) (or there were no more to be opened)
	 */
	bool Open();
	/**
	 * @brief Close the currently opened files
	 */
	void Close();

	/**
	 * @brief Used to reset the pointer to the first snp. If you need to scan through again, simply reset
	 */
	void Reset();

	/**
	 * @brief Returns the number of snps found in the buffer
	 */
	uint GetSnpCount();

	/**
	 * @brief returns the number of individuals in the file
	 */
	uint GetIndividualCount();

	/**
	 * @brief Returns the number of affected individuals have been identified
	 */
	//uint GetAffectedCount();
	
	/**
	 * @brief Returns the count of unaffected individuals identified
	 */
	//uint GetUnaffectedCount();

	 /**
	  * @brief Compiles the affected mask.	
 	  */ 
	void GetStatusMask(CaseControlStatus *);			
	void AppendScrambledStatus(CaseControlStatus *);

	/**
	 * Returns the line parser used to aquire the data
	 */
	GtLineParser *GetLineParser();
	/**
	 * @brief Produces the report for the effort performed
	 */
	void GenerateReport(ostream &os);

protected:
	StringArray filenames;						///<the array of input files
	GtLineParser *parser;						///<The line parser used buffer the input and produce snps
	uint currentSnp;							///<The snp index to be returned with NextSnp
};

inline
GtLineParser *GtFileParserBuffered::GetLineParser() {
	return parser;
}

inline
GtFileParserBuffered::GtFileParserBuffered(GtLineParser *parser, const char *filename) : parser(parser), currentSnp(0) {
	filenames.push_back(string(filename));
}

inline
GtFileParserBuffered::~GtFileParserBuffered()
{
	delete parser;
}



inline
void GtFileParserBuffered::GetStatusMask(CaseControlStatus *container) {
	parser->GetStatusMask(container);
}

inline
void GtFileParserBuffered::AppendScrambledStatus(CaseControlStatus *container) {
	parser->AppendScrambledStatus(container);
}

inline
uint GtFileParserBuffered::GetSnpCount() {
	return parser->GetSnpCount();
}

inline
uint GtFileParserBuffered::GetIndividualCount() {
	return parser->GetIndividualCount();
}
	


inline
SnpAligned *GtFileParserBuffered::NextSnp() {
	return parser->GetSnp(currentSnp++);
}



}

}

#endif
