//
// C++ Implementation: gtfileparserbuffered
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "gtfileparserbuffered.h"

namespace Genetics {

namespace Parser {

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

void GtFileParserBuffered::GenerateReport(ostream &os) {
	os<<setw(45)<<"Data Source: ";
	uint count=filenames.size();
	os<<"["<<filenames[0];
	for (uint i=1;i<count; i++) 
		os<<", "<<filenames[i];
	os<<"]\n";
	parser->GenerateReport( os );
}



void GtFileParserBuffered::Reset() {
	currentSnp=0;
}

void GtFileParserBuffered::Close() {
	Reset();
}


}

}
