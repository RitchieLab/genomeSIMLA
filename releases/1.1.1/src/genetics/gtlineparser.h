//
// C++ Interface: gtlineparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GENETICS_PARSERGTLINEPARSER_H
#define GENETICS_PARSERGTLINEPARSER_H
#include "snpaligned.h"
#include "snppool.h"

namespace Genetics {

namespace Parser {

typedef vector<vector<int > > MdrDataVector;

/**
	@brief Base class for the various line parsers. 
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class GtLineParser : public AsciiParser {
public:
	GtLineParser(int genoCount = 4);
	virtual ~GtLineParser();
	virtual SnpAligned *GetSnp(uint snpID)=0;			///<Returns the Snp at idx, snpID
	virtual uint GetIndividualCount()=0;				///<Returns the number of individuals found
	virtual uint GetSnpCount()=0;						///<Returns the number of snps encountered
	virtual void GenerateReport(ostream &os)=0; 		///<Generates a string report
	virtual AsciiParser *GetLineCounter();				///<Returns the object used for counting lines
	virtual void GetStatusMask(CaseControlStatus *)=0;			///<The container containing each of the Masks
	virtual void AppendScrambledStatus(CaseControlStatus *)=0;///<return a container for a single set of randomized tests
	virtual void PostParse();							///<Override this if the parser needs to do any work after the loading is completed
	
protected:
	BasicLineCounter lineCounter;
	int genoCount;										///<The number of genotypes expected

};


/** 
 * @brief Line Parser to be used with Base Pair style input files
 * @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class GtLineParserBP: public GtLineParser {
public:
	GtLineParserBP(uint affected, uint unaffected, int maxGenotypes);
	~GtLineParserBP();
	
	bool ParseLine(const char *line, uint val);		///<This is used to parser a single line from the line parser
	SnpAligned *GetSnp(uint snpID);						///<Returns the Snp at idx, snpID
	void GetStatusMask(CaseControlStatus *);				///<The bit Mask produced by the affected status
	void AppendScrambledStatus(CaseControlStatus *);		///<Returns the scrambled Mask for the affected status
	uint GetIndividualCount();							///<Returns the number of individuals found
	uint GetSnpCount();									///<Returns the number of snps encountered
	//uint GetAffectedCount();
	//uint GetUnaffectedCount();
	void GenerateReport(ostream &os);
protected:	
	CaseControlStatus statusMask;							///<Populated with the affected pattern
	vector<SnpAligned *> snps;							///<Storage for snps
	int individualCount;								///<the number of individuals encountered
	uint snpCount;										///<The number of snps encountered
	SnpPool *pool;										///<The pool used to acquire and release the snps
};



/** 
 * @brief Line Parser to be used with MDR formatted input files
 * @author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
 */
class GtLineParserMDR : public GtLineParser {
public:
	GtLineParserMDR(int affValue, bool readLabel=false, int genocount=4);
	~GtLineParserMDR();
	
	bool ParseLine(const char *line, uint val);		///<This is used to parser a single line from the line parser
	SnpAligned *GetSnp(uint snpID);						///<Returns the Snp at idx, snpID
	uint GetIndividualCount();							///<Returns the number of individuals found
	uint GetSnpCount();									///<Returns the number of snps encountered

	bool InitData(const char *line);					///<Initialize the snps based on a possible header

	int ConvertData(const char *value);					///<Used to convert the numerical values found in the mdr file
	void GetStatusMask(CaseControlStatus *);				///<The bit Mask produced by the affected status
	void AppendScrambledStatus(CaseControlStatus *);		///<Returns the scrambled Mask for the affected status
	void PrintData();									///<Print the contents of the array for debugging purposes
	void GenerateReport(ostream &os);
protected:
	/**
 	 * Status for this object is singular. We should only have a single entry for now.
	 */
	//StatusContainer *statusMask;						///<Populated with the affected pattern
	CaseControlStatus stMask;
	bool readLabel;										///<Indicates if there is expected to be a header
	MdrDataVector data;									///<The array of integers
	vector<SnpAligned::Label> labels;					///<The array of labels
	int affectedValue;									///<If our data has unknown status, this has to be changed - This seems to be acceptable, but users should be aware to clean their data before hand
	int unknownGenotype;								///<This is the value used to assign
	SnpPool *pool;										///<Snp pool used to acquire and release items
};



//---------------------------GtLineParser

inline
GtLineParser::GtLineParser(int genoCount) : genoCount(genoCount) {}

inline
GtLineParser::~GtLineParser() { }

inline
void GtLineParser::PostParse() { }

inline
AsciiParser *GtLineParser::GetLineCounter() {
	return &lineCounter;
}

//---------------------------GtLineParserBP 

inline
GtLineParserBP::GtLineParserBP(uint affected, uint unaffected, int maxGenotypes) : GtLineParser(maxGenotypes), statusMask(0), individualCount(-1), snpCount(0)  {
	pool=SnpPool::Instance();
	BitSetType affectedMask(affected+unaffected, false);
	
	for (uint i=0; i<affected; i++) 
		affectedMask[i]=true;
	BitSetType unaffMask = ~affectedMask;
	statusMask=CaseControlStatus(affectedMask, unaffMask);
 }


inline
GtLineParserBP::~GtLineParserBP() { 
	vector<SnpAligned *>::iterator i=snps.begin();
	vector<SnpAligned *>::iterator end=snps.end();

	SnpAligned *snp;
	for (; i!=end; ++i) {
		snp=*i;
		pool->ReleaseSnp(snp);
	}
}



inline
uint GtLineParserBP::GetIndividualCount() {
	return individualCount;
}

inline
uint GtLineParserBP::GetSnpCount() {
	return snpCount;
}

inline
SnpAligned *GtLineParserBP::GetSnp(uint snpID) {
	SnpAligned *snp=NULL;
	
	if (snpID<snps.size()) {
		snp=snps[snpID];
		snp->IncrementInstanceCount();			///<If you request a snp more than once, you better release it that many times!
	}
	return snp;
}

inline
bool GtLineParserBP::ParseLine(const char *line, uint val) {
	SnpAligned *snp=pool->GetSnp(genoCount);
	int indCount = snp->ImportBpSnp( val, ++snpCount, line);

	if (indCount>0)	{
		if (individualCount == -1)
			individualCount=indCount;
		else	
			if (individualCount!=indCount) {
				cout<<"!! Unbalanced data encountered on line "<<val<<". Previously, "<<individualCount<<" but "<<indCount<<" were found. ("<<line<<")\n";
				assert(individualCount==indCount);
			}
		snps.push_back(snp);
		return true;
	}	
	else	{
		pool->ReleaseSnp(snp);
		return false;
	}
}


inline
uint GtLineParserMDR::GetSnpCount() {
	if (data.size() > 0) 
		return data[0].size() - 1;
	else
		return 0;
}

inline
uint GtLineParserMDR::GetIndividualCount() {
	return data.size();
}




inline
GtLineParserMDR::GtLineParserMDR(int affValue, bool readLabel, int genotypeCount/*=-1*/) : GtLineParser(genotypeCount), stMask(0), readLabel(readLabel), affectedValue(affValue), unknownGenotype(-1) {
	pool=SnpPool::Instance();
}
inline
GtLineParserMDR::~GtLineParserMDR() {
}



}
}
#endif
