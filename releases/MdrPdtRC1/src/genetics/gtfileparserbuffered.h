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
void GtFileParserBuffered::GenerateReport(ostream &os) {
	os<<setw(45)<<"Data Source: ";
	uint count=filenames.size();
	os<<"["<<filenames[0];
	for (uint i=1;i<count; i++) 
		os<<", "<<filenames[i];
	os<<"]\n";
	parser->GenerateReport( os );
}

/**
inline
uint GtFileParserBuffered::GetAffectedCount() {
	return parser->GetAffectedCount();
}

inline
uint GtFileParserBuffered::GetUnaffectedCount() {
	return parser->GetUnaffectedCount();
}
*/
inline
void GtFileParserBuffered::GetStatusMask(CaseControlStatus *container) {
	parser->GetStatusMask(container);
}

inline
void GtFileParserBuffered::AppendScrambledStatus(CaseControlStatus *container) {
	parser->AppendScrambledStatus(container);
}

inline
void GtFileParserBuffered::Reset() {
	currentSnp=0;
}
inline
void GtFileParserBuffered::Close() {
	Reset();
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
bool GtFileParserBuffered::Open() {
	uint count = 0;		
	int fileCount = filenames.size();
	string filename;
	//AsciiParser *lineCounter=parser->GetLineCounter();
	
	LineParser lp;
	//For now, let's just look at the first file in any multi file bunch
	for (int i =0; i<fileCount; i++) {
		filename=filenames[i];
//		count += lp.Parse( filename.c_str(), lineCounter);
		count += lp.Parse( filename.c_str(), parser);
	}
	
	parser->PostParse();

	return count;	
}

inline
SnpAligned *GtFileParserBuffered::NextSnp() {
	return parser->GetSnp(currentSnp++);
}



}

}

#endif
