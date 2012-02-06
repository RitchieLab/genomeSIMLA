//
// C++ Implementation: ldplotmanager.h
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldplotmanager.h"
#include "ldprimepngwriter.h"
#include <iomanip>
#include <boost/timer.hpp>
#include <sstream>
#include "ldlineplot.h"


namespace Simulation {
namespace Visualization {



//We are assuming that the allele frequencies are already calculated prior to establishing the plotter object
LdPlotManager::LdPlotManager(vector<Locus> * loci) : results(NULL), loci(loci) { }


LdPlotManager::~LdPlotManager() {


}

void LdPlotManager::Init(vector<vector<LdResult> >& results, BlockListHead<BlockListNodeFourGammetes>* haplotypeBlocks) {
	this->results = &results;
	haplotypes = haplotypeBlocks;
}


void LdPlotManager::SetLabel(const char *lbl) {
	label = lbl;
}

string LdPlotManager::WriteLdDPrime(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl) {
	char adjustedFilename[1024];
	sprintf(adjustedFilename, "%s%s-dp", filename, lbl);

	//The results contain loci IN ADDITION to the snp of interest
	uint lociCount = results->size() + 1;
	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	size_t first = 0;
	if (block)
		first = block->GetFirstIndex();
	

	if (first > ldBufferSize)
		first-=ldBufferSize;
	else
		first = 0;

	size_t last = lociCount - 1;
	if (block)
		last = block->GetLastIndex() + ldBufferSize;

	if (last >= lociCount)
		last = lociCount - 1;

	assert(last<lociCount);

	//Since we want to limit the size of the plot to heights that will actually be used, we need to 
	//scan over all possible couplings to determine who is within range for LD consideration
	size_t cols, rows;
	int firstLocation, lastLocation;
	CountRC(first, last, rows, cols, firstLocation, lastLocation);
	

	string pngFilename = "";

	if (rows && cols) {
		//We'll set up the image parameters for the current region.
		ImageParameters *params = new ImageParameters;
		params->Init(cols, rows);

		LdWriter *png = new LdPrimePngWriter(params, firstLocation, lastLocation);

		stringstream lbl;
		lbl<<label;
		if (block)
			lbl<<"("<<block->GetLabel()<<")";
		lbl<<": D'";
		png->SetLabel(lbl.str().c_str());

		//((LdDPrimePng *)png)->SetLabel(label.c_str());
		imgHeight = params->dimensions.y;
		imgWidth  = params->dimensions.x;

		//We'll write the LD values to file
		pngFilename = WriteLdData( adjustedFilename, first, last, png, params->snpDepth);


		//Then add each of the blocks- making the central snp more boldly marked
		blocks = haplotypes->GetValidBlocks();
		haplotypes->GetBlockDensities(minBlockDensity, maxBlockDensity);

		((LdPrimePngWriter *)png)->SetBlockDensity(minBlockDensity, maxBlockDensity);
		BlockList::iterator itr=blocks->begin();
		BlockList::iterator end=blocks->end();

		while (itr != end) {
			png->AddBlock(*(itr++), 1);
		}
		if (block)
			png->AddBlock( block, 3 );
		png->Close();
		delete png;
	}
	else {
		cout<<"Unable to write an empty png!\n";
		abort();
	}
	return pngFilename;
}
//string LdPlotManager::WriteLdRSquared(const char *filename, HaplotypeBlock *block, int &imgHeight, int &imgWidth, uint ldBufferSize, const char *lbl) {
string LdPlotManager::WriteLdRSquared(const char *filename, BlockListNode *block, uint &imgHeight, uint &imgWidth, uint ldBufferSize, const char *lbl) {
	char adjustedFilename[1024];
	sprintf(adjustedFilename, "%s%s-rs", filename, lbl);

	uint lociCount = results->size();
	assert(results->at(lociCount-1).size() == 1);

	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	size_t first = 0;
	if (block)
		first = block->GetFirstIndex();

	if (first > ldBufferSize)
		first-=ldBufferSize;
	else
		first = 0;

	size_t last = lociCount - 1;
	if (block)
		last = block->GetLastIndex() + ldBufferSize;

	if (last >= lociCount)
		last = lociCount - 1;

	assert(last<lociCount);

	//Calculate the number of snps to be considered as well as the largest pairing based on the max distance setting
	size_t cols, rows;
	int firstLocation, lastLocation;
	CountRC(first, last, rows, cols, firstLocation, lastLocation);

	string pngFilename = "";

	if (rows && cols) {
		//Create the image parameters and set up the renderer 
		ImageParameters *params = new ImageParameters;
		params->Init(cols, rows);
		LdWriter *png = new LdRSquaredPngWriter(params, firstLocation, lastLocation);
		((LdRSquaredPngWriter *)png)->SetBlockDensity(minBlockDensity, maxBlockDensity);
		//((LdDPrimePng *)png)->SetLabel(label.c_str());
		imgHeight = params->dimensions.y;
		imgWidth  = params->dimensions.x;

		stringstream lbl;
		lbl<<label;
		if (block)
			lbl<<"("<<block->GetLabel()<<")";
		lbl<<": R^2";
		png->SetLabel(lbl.str().c_str());

		pngFilename = WriteLdData( adjustedFilename, first, last, png, params->snpDepth);

		//Iterate over all blocks and add them to the plot. If it's the central block, we'll make it darker
		BlockList::iterator itr=blocks->begin();
		BlockList::iterator end=blocks->end();
		while (itr != end) {
			png->AddBlock(*(itr++), 1);
		}
		if (block)
			png->AddBlock( block, 3 );
		png->Close();
		delete png;
	}
	else {
		cout<<"Unable to write an empty png!\n";
		abort();
	}
	return pngFilename;
}

string LdPlotManager::WriteLdReport(const char *filename, BlockListNode*block) {
	LdWriter *txtReport = new LdTextReport();
	string r = WriteLdData(filename, block, txtReport, 0);
	delete txtReport;
	return r;
}


string LdPlotManager::WriteLdData(const char *filename, BlockListNode *block, LdWriter* writer, uint padding) {
	uint lociCount = results->size();

	//Make sure we have reasonable limits. The last value should be inclusive (so if it is 9, loci[9] must be meaningful
	int first = 0;
	int last = lociCount - 1;
	if (block) {
		first = block->GetFirstIndex() - padding;
		last = block->GetLastIndex() + padding;// + ldBufferSize;
	
		if (first < 0)
			first = 0;
		if (last >= lociCount)
			last = lociCount - 1;
	}

	assert(first<lociCount);
	assert(last<lociCount);

	return  WriteLdData(filename, first, last, writer, 0);		
}
 
/**
 * @NOTE it is very important to realize that the cache is not smart at all. If you initialize the LD Cache to a 
 *       subset of the Loci, and later attempt to render the whole chromosome, you will likely get a crash. Definitely, you
 *       will get questionable results. 
 */
string LdPlotManager::WriteLdData(const char *filename, uint first, uint last, LdWriter* writer, uint maxSnpDepth) {
	//This will open the "writer" object and write up any header information
	string pngFilename = writer->Open(filename, *loci, first, last);

	assert(results);

	int snpCount = results->size();
	vector<vector<LdResult> >::iterator itr = results->begin();
	vector<vector<LdResult> >::iterator end = results->end();

	//Move up to the first SNP of interest
	int idx = 0;
	while (idx++ < first)
		itr++;
	idx = first;
	
	while (idx < last && itr != end) {
		vector<LdResult>& resultArray = *itr;
		vector<LdResult>::iterator resItr = resultArray.begin();
		vector<LdResult>::iterator resEnd = resultArray.end();
		Locus &l = loci->at(idx);
		writer->WriteHeader(idx - first, l);
//cerr<<"Locus("<<idx-first<<") "<<l<<"[";

		int startIdx = idx;
		for (; resItr != resEnd && startIdx++ < last; resItr++ ) {
			LdResult &result = *resItr;
//cerr<<" "<<result.Q->GetLabel();
			writer->Write(*(result.Q), result.dprime, result.lod, result.rsquared);
		}
//cerr<<" ]\n";
		idx++;
		itr++;
	}

	return pngFilename;
}

void LdPlotManager::CountRC(size_t first, size_t last, size_t &rows, size_t &cols, int &firstLocation, int &lastLocation) {
	cols = last-first+1;
	//Grab the appropriate set of results from the overall container
	vector<LdResult>& resultArray = results->at(first);
	rows=resultArray.size();
	if (rows > 0) {
		//Get the location of the first snp
		Locus &l1 = loci->at(first);
//Locus *l1 = resultArray[0].P;
		firstLocation = l1.GetLocation();	//resultArray[0].P->GetLocation();
assert(firstLocation == resultArray[0].P->GetLocation());
		vector<LdResult>& endArray = results->at(first + cols-2);
//		assert(endArray.size() > 0);
		Locus &l2 = loci->at(last);
//Locus *l2 = endArray[0].Q;
		lastLocation = l2.GetLocation();
if (endArray.size() > 0) {
assert(lastLocation == endArray[0].Q->GetLocation());
}
//		lastLocation=endArray[0].Q->GetLocation();
	}
	vector<vector<LdResult> >::iterator itr = results->begin();
	vector<vector<LdResult> >::iterator end = results->end();

	int idx = 0;
	while (itr != end) {
		if (idx >= first && idx <= last){ 
			if (itr->size() > rows)
				rows = itr->size();
		}
		idx++;
		itr++;
	}	
	if (rows > cols)
		rows = cols;
}

}
}
