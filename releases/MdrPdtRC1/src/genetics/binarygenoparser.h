//
// C++ Interface: binarygenoparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ESEBINARYGENOPARSER_H
#define ESEBINARYGENOPARSER_H
#include <fstream>
#include "snppool.h"
#include "snpaligned.h"

namespace Genetics {

/**
@brief Handles the parsing of genotype information in a compressed binary representation. 
<B>Usage:</B> Upon opening, the header data will be read and the details about the file can be accessed
at any time using the appropriate accessors. Once the file has successfully been opened, calls can be
made to GetNextSnp() which returns the snp found next in the file or NULL once the end of file has been
encountered.

<P><B>File Format:</B>
<UL>
<LI>File description containing versioning details for the application that produced the file
<LI>File Type (unsigned long)
<LI>Genotype Count (int)
<LI>Case Count (unsigned long)
<LI>Control Count (unsigned long)
<LI>Label Style (usigned short)
</UL>
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class BinaryGenoParser{
public:
    BinaryGenoParser(const char *filename);

    ~BinaryGenoParser();

	/**
	 * @brief Opens the file and parses the header for details about the contents.
	 */
	bool Open();
	
	/**
	 * @brief Parses next snp in the file.
	 * @return Valid Snp or NULL
	 */
	SnpAligned *GetNextSnp();

	/**
	 * @brief Set the width of the actual encoded value.
	 * By default, the value is set to 2 bits. If there are more than 4 genotypes, this will need to 
	 * be increased to accomodate the number of genotypes. 
	 * This should probably become automatic
	 */
	void SetFramesize(uint size);
	void SetWordSize(uint size);

	string GetFilename();
	uint GetFileType();
	uint GetSnpCount();
	uint GetGenoCount();
	uint GetCaseCount();
	uint GetControlCount();
	uint GetLabelStyle();
	uint GetBlockSize();
	

protected:
	ifstream dataStream;					///<This is the input stream containing the snps
	string filename;						///<The file that is associated with this object
	uint fileType;							///<Figure out the values and meaning of this
	uint snpCount;							///<The number of snps encountered thus far
	uint genoCount;							///<The number of genotypes found so far
	uint caseCount;							///<The number of cases
	uint controlCount;						///<The number of controls
	uint labelStyle;						///<I am not sure what this does yet! It's in the file header, though
	SnpPool *pool;							///<Used to acquire the snps to be used

	uint blockSize;							///<The number of storage units required to save the bits
	uint frameSize;							///<This is how many bits is used to store the integer information in the compressed string
	uint wordSize;							///<The size of a single unit of storage. Default value, 32
	uint curSnp;							///<This is used to help make sure we don't read too many

};

inline
BinaryGenoParser::BinaryGenoParser(const char *filename) : 
	filename(filename),
	fileType(0),
	snpCount(0),
	genoCount(4),
	caseCount(0),
	controlCount(0),
	labelStyle(0),
	frameSize(2),
	wordSize(32),
	curSnp(0)
{
	pool = SnpPool::Instance();
		
}	

inline
void BinaryGenoParser::SetWordSize(uint size) {
	wordSize=size;
}

inline
void BinaryGenoParser::SetFramesize(uint size) {
	frameSize=size;
	genoCount=2^size;
}

inline
uint BinaryGenoParser::GetBlockSize() {
	return blockSize;
}

inline
string BinaryGenoParser::GetFilename() {
	return filename;
}

inline
uint BinaryGenoParser::GetFileType() {
	return fileType;
}

inline
uint BinaryGenoParser::GetSnpCount() {
	return snpCount;
}

inline
uint BinaryGenoParser::GetCaseCount() {
	return caseCount;
}

inline
uint BinaryGenoParser::GetControlCount() {
	return controlCount;
}

inline
uint BinaryGenoParser::GetLabelStyle() {
	return labelStyle;
}



}

#endif
