//
// C++ Implementation: ldwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldwriter.h"
#include <iomanip>
#include <math.h>

namespace Simulation {

namespace Visualization {


LdWriter::LdWriter() : A(NULL), curSnpIdx(0)
{
}


LdWriter::~LdWriter()
{
}


string LdTextReport::Open(const char *filename, vector<Locus> &loci, uint first, uint last) {
	char ldFilename[1024];
	char mafFilename[1024];

	sprintf(ldFilename, "%s.ld", filename);
	sprintf(mafFilename, "%s.maf", filename);

	ldFile.open(ldFilename, ios::out);
	maf.open(mafFilename, ios::out);

	//First column is SNP. Second is MAF, ranged 0.0-0.5
	maf<<"SNP\tMAF\t0.0\t0.5\n";
	ldFile<<"L1\tL2\tD'\tLOD\tr2";
	ldFile<<"\n";
	return ldFilename;
}
void LdTextReport::Close() {
	ldFile.close();
	ldFile.close();
}
void LdTextReport::WriteHeader(uint idx, Locus &A) {
//	maf<<A.GetLabel()<<"\t"<<A.GetMinAlleleFreq()<<"\n";
	this->A = &A;
	curSnpIdx++;	
}
void LdTextReport::Write(Locus &B, float &dprime, float &lod, float &rsquared) {
	ldFile<<A->GetLabel()
		<<"\t"<<B.GetLabel()
		<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<dprime
		<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<lod
		<<"\t"<< setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<rsquared
		<<"\n";

}



}
}
