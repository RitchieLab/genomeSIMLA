//
// C++ Implementation: binarygenoparser
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "binarygenoparser.h"
#include "utility/utility.h"
namespace Genetics {

BinaryGenoParser::~BinaryGenoParser()
{
}


/**
 * @brief Opens the file and parses the header for details about the contents
 */
bool BinaryGenoParser::Open() {
	dataStream.open(filename.c_str(), ios_base::in);

	if (dataStream.eof() || dataStream.fail() || dataStream.bad() ) {
		return false;
	}
	
	//First, let's skip over the comment at the beginning of the file
	//dataStream.read(line, 1000);


	//Read in header details
	dataStream.read((char*)&fileType, sizeof(unsigned long int));
	dataStream.read((char*)&snpCount, sizeof(int));
	dataStream.read((char*)&caseCount, sizeof(unsigned long int));
	dataStream.read((char*)&controlCount, sizeof(unsigned long int));
	dataStream.read((char*)&labelStyle, sizeof(unsigned short));

	curSnp=0;

	/**
 	 * @todo This should be something flexible, since the actual bin array parser is 
 	 */
	//Figure out how many blocks we need to read in
	uint indivCount=caseCount+controlCount;
	if (indivCount%wordSize == 0) 
		blockSize=indivCount/wordSize;
	else
		blockSize=indivCount/wordSize+1;



	return true;
}
	
/**
 * @brief Parses next snp in the file.
 * @return Valid Snp or NULL
 */
SnpAligned *BinaryGenoParser::GetNextSnp() {
	SnpAligned *snp=NULL;
	uint indCount = caseCount+controlCount;

	if (curSnp < snpCount && !dataStream.eof()) {
		uint label;
		uint locusToLineNumber;
		unsigned char genotypeValues;
		//uint labelStyle;
		char trash[4];
	
		//Read in the snp descriptor details
		dataStream.read((char*)&label, sizeof(uint));
		dataStream.read((char*)&locusToLineNumber, sizeof(uint));
		dataStream.read((char*)&genotypeValues, sizeof(unsigned char));
		dataStream.read((char*)&trash, 3);									///<Alignment problem

		//If this get's tripped, then there is a difference between the two. I don't like carrying around this
	 	//exta integer with these little snps. We can keep up with ranges from files through some other means
		if (label!=locusToLineNumber)
			cout<<"-- Label and line number don't match! "<<label<<" - "<<locusToLineNumber<<"\n";
		assert(label == locusToLineNumber);
	
		uint *data = new uint[blockSize];
		uint a;
		//Loads the entire block
		//dataStream.read((char*)&data, sizeof(uint)*blockSize);
		for (uint i=0; i<blockSize; i++) {
			dataStream.read((char*)&a, sizeof(uint));

			data[i]=a;
		}
		snp=pool->GetSnp(genoCount);
		snp->SetIndividualCount(indCount);
		//cout<<"Label: "<<label<<"\tLocus: "<<locusToLineNumber<<"\tGenotypes: "<<genotypeValues<<"\t";
		curSnp++;
		snp->ImportCompressedBinSnp(locusToLineNumber, label, data, blockSize, frameSize);
		cout<<"Snp "<<label<<" : "<<snp->toString()<<"\n";
		delete[] data;
	}		
	return snp;
}


}
