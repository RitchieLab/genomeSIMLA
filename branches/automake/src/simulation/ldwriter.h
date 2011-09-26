//
// C++ Interface: ldwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDWRITER_H
#define SIMULATION_VISUALIZATIONLDWRITER_H
#include <fstream>
#include <string>
#include "locus.h"
#include "pngwriter/pngwriter.h"
//#include "haplotypeblock.h"
#include "blocklistnodefourgammetes.h"

namespace Simulation {

namespace Visualization {

using namespace std;



/**
@brief various classes associated with writing LD Data

	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdWriter{
public:
    LdWriter();
    virtual ~LdWriter();
	virtual void WriteHeader(uint idx, Locus &A) = 0;
	virtual void Write(Locus &B, float &dprime, float &lod, float &rsquared) = 0;
	virtual string Open(const char *filename, vector<Locus> &loci, uint first, uint last) = 0;
	virtual void Close() =0;
	virtual void AddBlock(BlockListNode *block, uint classification) = 0;
	virtual void SetLabel(const char *) { cout << "ASDFASDFASDFASDF";}
//AddBlock(HaplotypeBlock *block, uint classification)=0;
protected:
	Locus *A;
	int curSnpIdx;
};

class LdTextReport : public LdWriter {
	string Open(const char *filename, vector<Locus> &loci, uint first, uint last);
	void WriteHeader(uint idx, Locus &A);
	void Write(Locus &B, float &dprime, float &lod, float &rsquared);
	void Close();
	void AddBlock(BlockListNode *block, uint classification) {}
//AddBlock(HaplotypeBlock *block, uint classification) {};
protected:
	ofstream ldFile;
	ofstream maf;
	string lastLocus;
};

class LdWriterPng : public LdWriter {
public:
//	static int spacerSize;
	static int standardBlockSize;
	static int tinyBlockSize;
	static int medBlockSize;
	static int marginX;
	static int marginY;
	static string font;
};




}

}

#endif
