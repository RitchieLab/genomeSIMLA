//
// C++ Implementation: blocklistnodefourgammetes
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "blocklistnodefourgammetes.h"
#include <iomanip>
#include "utility/strings.h"

namespace Simulation {

using namespace Utility;

BlockListNodeFourGammetes::BlockListNodeFourGammetes(float criterion, BlockListNode *previous, BlockList *completedList, BlockList *deletedList) :
		BlockListNode(criterion, previous, completedList, deletedList), rowsWritten(0), maxImageWidth(1024) {

}

BlockListNodeFourGammetes::BlockListNodeFourGammetes(float criterion) : BlockListNode(criterion) {
}


BlockListNodeFourGammetes::~BlockListNodeFourGammetes()
{
//	cout<<"Deleting Block: "<<firstIdx<<" : "<<lastIdx<<"\n";
}


BlockListNode *BlockListNodeFourGammetes::NewNode(BlockListNode *previous, uint firstIdx, uint lastIdx) {
	BlockListNode *newNode;
	if (deletedList->size() > 0) {
		newNode = deletedList->back();
		deletedList->pop_back();
		((BlockListNodeFourGammetes*)newNode)->previous = previous;
	}
	else {
		newNode = new BlockListNodeFourGammetes(criterion, previous, completedList, deletedList);
	}
	newNode->SetIndexBounds(firstIdx, lastIdx);

	return newNode;
}

bool BlockListNodeFourGammetes::ValidateBlocks(vector<Locus>& sourceLoci, BlockFilter *filter) {
	//Calculate avg MAF
	float sumMAF = 0.0;
	firstLocation = sourceLoci[firstIdx].GetLocation();
	lastLocation  = sourceLoci[lastIdx].GetLocation();
	totalDistance = lastLocation -  firstLocation;
	loci.clear();

	if (filter)
		filter->Reset();

	float minMAF = 0.6;
	float maxMAF = 0.0;
	for (size_t cur = firstIdx; cur<=lastIdx; cur++) {
		Locus &l = sourceLoci[cur];
		float maf = l.GetMinAlleleFreq();
		sumMAF+=maf;
		if (maf < minMAF) {
			minMAF = maf;
			minAlleleFreq = &l;
		}
		if (maf > maxMAF) {
			maxMAF = maf;
			maxAlleleFreq = &l;
		}
		loci.push_back(l);
		if (filter)
			filter->Evaluate(l);
	}
	avgMAF = sumMAF / (float)loci.size();
	blockDensity = totalDistance / loci.size();
	if (filter) {
		selectionScore = filter->GetScore();
		return filter->IsValid();
	}
	//If we don't have a filter, we'll accept anything and they'll all have the same score
	return true;
}

/**
 * @brief generate summary report
 */
void BlockListNodeFourGammetes::Report(ostream& os) {
	size_t end=loci.size() - 1;
	os<<"Block: "<<loci[0].GetLabel()<<" - "<<loci[end].GetLabel()<<"  Size: "<<end+1<<" Average MAF: "<<GetAverageMAF()<<" BD: "<<GetBlockDensity()<<"\n";
}

/**
 * @brief Generate a more detailed report of the block contents
 */
void BlockListNodeFourGammetes::DetailedReport(ostream &os) {
	//OK, for ease of reading, we'll sort them- but we want to make a copy first
	size_t lociCount = this->loci.size();
	vector<Locus>loci = this->loci;

	//sort(loci.begin(), loci.end(), LTMAF());
	os<<"\tHaplotype Block, Size: "<<lociCount<<" - "<<GetBlockDensity()<<"\n";
	os<<"\t\tBlock Constituents: "
			<<minAlleleFreq->GetLabel()<<" - "
			<<maxAlleleFreq->GetLabel()<<" Avg MAF: "
			<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<GetAverageMAF()<<"\n";
    os<<setw(12)<<right<<"Index"<<setw(12)<<right<<"Label"<<setw(10)<<"Location"<<setw(12)<<setprecision(5)<<"Min. Al. Freq"<<"\n";
	for (size_t i=0; i<lociCount; i++) {
		Locus &loc = loci[i];
    	os<<setw(12)<<right<<i + 1<<
			setw(12)<<right<<loc.GetLabel()<<
			setw(10)<<loc.GetLocation()<<
			setw(12)<<setprecision(5)<<loc.GetMinAlleleFreq()<<"\n";	
	}
	
}


	
	
/**
 * @brief This is used for sorting in STD sort functions
 * @note We want it to work out so the best block is the first in the vector
 */
//bool BlockListNodeFourGammetes::operator<(const HaplotypeBlock& other) const {
bool BlockListNodeFourGammetes::operator<(const BlockListNodeFourGammetes& other) const {
	if (selectionScore == other.selectionScore)
		return loci.size() < other.loci.size();
	else
		return selectionScore < other.selectionScore;
}





static string tableHeaderColor = "#000000";
void BlockListNodeFourGammetes::HtmlSummaryBegin(ostream &os) {
	os<<"\n<DIV id=\"HaplotypeSummaryTable\"><!Summary Table for remaining blocks>\n";
	//os<<"<TABLE cellspacing=\"1\" cellpadding=\"1\" border=\"1\"><CAPTION>Block Summary</CAPTION>\n";
	os<<"<TABLE><CAPTION>Block Summary</CAPTION>\n";
	os<<"\t<tr><TH>Bounding SNPs</TH><th>Block Size</th><th>Avg. Min. Allele Freq.</TH><TH>Block Density</TH></STRONG></tr>\n";
}
void BlockListNodeFourGammetes::HtmlSummary(ostream& os, const char *linkToDetails) {
	size_t end=loci.size() - 1;

	static string rowBackgrounds[] = {"\t<TR><TD>", "\t<TR class=\"invert\"><TD>"};

	int idx = rowsWritten++ %2;
	os<<rowBackgrounds[idx];
	if (linkToDetails) 
		os<<"<A HREF=#"<<linkToDetails<<">"<<loci[0].GetLabel()<<"_"<<loci[end].GetLabel()<<"</A>";
	else
		os<<"<A name='"<<loci[0].GetLabel()<<"_"<<loci[end].GetLabel()<<"'>"
		  <<loci[0].GetLabel()<<"_"<<loci[end].GetLabel()<<"</A>";
	os<<"</TD><TD>"<<end+1
		<<"</TD><TD>"<<GetAverageMAF()
		<<"</TD><TD>"<<GetBlockDensity()
		<<"</TD><TR>\n";
}

void BlockListNodeFourGammetes::HtmlSummaryEnd(ostream &os) {
	os<<"</TABLE></DIV>\n\n";
}

string BlockListNodeFourGammetes::HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails) {
	size_t lociCount = this->loci.size();
	vector<Locus>loci = this->loci;
	string anchor=loci[0].GetLabel()+"_"+loci[lociCount - 1].GetLabel();
	//Let's get rid of the "/"s in the filename

	size_t width = maxImageWidth;
	if (imgWidth < maxImageWidth)
		width = imgWidth;

	//size_t offset = GetFirstIndex();

	sort(loci.begin(), loci.end());		//, LTMAF());

	os<<"\n</DIV><DIV class='spacer'></DIV><DIV id=\"HaplotypeBlock"<<loci[0].GetLabel()<<"-"<<loci[loci.size() - 1].GetLabel()<<"\" width="<<width<<" class='image-frame'>\n<CENTER>";
	os<<"<H2><A name="<<anchor<<">Haplotype Block "<<loci[0].GetLabel()<<" to "<<loci[loci.size() - 1].GetLabel()<<"</A></H2>\n";

	if (strlen(dPrimeFilename) > 0) {
		const char *picFirst = strrchr(dPrimeDetails, '/');
		if (picFirst)
			picFirst++;
		else 
			picFirst = dPrimeDetails;
	
		const char *zoomedOut = strrchr(dPrimeFilename, '/');
		if (zoomedOut)
			zoomedOut++;
		else
			zoomedOut = dPrimeFilename;
		

		os<<"\t<IMG SRC='"<<ExtractFilename(picFirst)<<"' border=\"0\" alt=\"Haplotype Plot\" id=\""<<ExtractFilename(dPrimeFilename)<<"\" width="<<width<<" onClick=\"window.open('"<<ExtractFilename(zoomedOut)<<"')\"> </IMG>\n";
	}

	if (strlen(rSquaredFilename) > 0) {
		const char *picFirst = strrchr(rSquaredDetails, '/');
		if (picFirst)
			picFirst++;
		else
			picFirst = rSquaredDetails;
			
		const char *zoomedOut = strrchr(rSquaredFilename, '/');
		if (zoomedOut)
			zoomedOut++;
		else 
			zoomedOut = rSquaredFilename;
		os<<"\t<IMG SRC='"<<ExtractFilename(picFirst)<<"' border=\"0\" alt=\"Haplotype Plot\" id=\""<<ExtractFilename(rSquaredFilename)<<"\" width="<<width<<" onClick=\"window.open('"<<ExtractFilename(zoomedOut)<<"')\"> </IMG>\n";
	}		
	

	os<<"\t<TABLE><CAPTION>Block Details</CAPTION><TR><TH>Size</TH><TH>Block Density</TH><TH>Avg. MAF</TH></TR>\n";
	os<<"\t<TR class='invert'><TD>"<<lociCount<<"</TD> <TD>"<<GetBlockDensity()<<"</TD> <TD>"<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<GetAverageMAF()<<"</TD></TR>\n\t</TABLE>";

	os<<"\t<DIV class='spacer'></DIV><TABLE class=\"SummaryTable\"><CAPTION>Block Constituents</CAPTION>\n";
	os<<"\t\t<TR><TH>Index</TH><TH>Label</th><th>Location</th><th>Freq (Al 1)</TH><TH>Freq (Al 2)</TH></tr>\n";

	for (size_t i=0; i<lociCount; i++) {
		Locus &loc = loci[i];

		int idx = i%2;
		if (idx == 0)
			os<<"<TR class=\"invert\">";
		else
			os<<"<TR>";
	
    	os<<"<TD>"<<setw(12)<<right<<i<<"("<<loc.GetChromID() + 1<<":"<<loc.GetID() + 1<<")"<<
			"</TD><TD>"<<setw(12)<<right<<loc.GetLabel()<<
			"</TD><TD>"<<setw(10)<<loc.GetLocation()<<
			"</TD><TD>"<<setw(12)<<setprecision(5)<<loc.Freq1()<<"</TD>"<<
			"</TD><TD>"<<setw(12)<<setprecision(5)<<loc.Freq2()<<"</TD></TR>\n";	
	}
	os<<"\t</TABLE></DIV>\n";
	return anchor;
}



}
