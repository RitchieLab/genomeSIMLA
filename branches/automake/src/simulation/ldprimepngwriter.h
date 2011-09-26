//
// C++ Interface: ldprimepngwriter
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SIMULATION_VISUALIZATIONLDPRIMEPNGWRITER_H
#define SIMULATION_VISUALIZATIONLDPRIMEPNGWRITER_H

#include "ldpngcomponent.h"
#include "ldlabelregion.h"
#include "ldlocationbar.h"
#include "ldwriter.h"
#include "ldmafregion.h"
#include "ldfooterregion.h"

namespace Simulation {

namespace Visualization {

/**
	@author Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>
*/
class LdPrimePngWriter : public LdWriter {
public:
    LdPrimePngWriter(ImageParameters *params, uint start, uint stop);
    ~LdPrimePngWriter();
	
	string Open(const char *filename, vector<Locus> &loci, uint first, uint last);
	void Close();
	void Init(vector<Locus> &loci, uint first, uint last);
	virtual void Write(Locus &B, float &dprime, float &lod, float &rsquared);
//	void AddBlock(HaplotypeBlock *block, uint classification);
	void AddBlock(BlockListNode *block, uint classification);
	void WriteHeader(uint idx, Locus &A);
	void SetBlockDensity(uint min, uint max);
	ImageParameters *GetImageParameters() {
		return imgParams;
	}
	void SetLabel(const char *label) { footer.SetLabel(label); }

protected:
	ImageParameters *imgParams;
	LdLabelRegion labels;
	LdLocationBar locations;
	LdMAFRegion	  maf;
	LdFooterRegion footer;
	size_t idx1;
	size_t idx2;
	
	size_t firstSnpIndex;				///<First Snp of the rendered area
	size_t lastSnpIndex;				///<Last snp of the rendered area
	float bRed;							///<Block marker red
	float bGreen;						///<Block marker green
	float bBlue;						///<Block marker blue

};

class LdRSquaredPngWriter : public LdPrimePngWriter {
public:
	LdRSquaredPngWriter(ImageParameters *params, uint start, uint stop);
	~LdRSquaredPngWriter();
	virtual void Write(Locus &B, float &dprime, float &lod, float &rsquared);
};


}

}

#endif
