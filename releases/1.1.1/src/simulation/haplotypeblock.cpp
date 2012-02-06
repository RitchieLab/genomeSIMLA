//
// C++ Implementation: haplotypeblock
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
//#include "haplotypeblock.h"
#include <iomanip>

namespace Simulation {

uint HaplotypeBlock::rowsWritten = 0;
uint HaplotypeBlock::maxImageWidth = 1024;

HaplotypeBlock::HaplotypeBlock(size_t idx1, Locus *loc1, size_t idx2, Locus *loc2) : firstLocusIndex(-1), lastLocusIndex(0), minAlleleFreq(NULL), maxAlleleFreq(NULL), totalMAF(0.0), lastIndex(0) {

	Append(idx1, loc1);
	Append(idx2, loc2);	
}


HaplotypeBlock::~HaplotypeBlock() {

}

void HaplotypeBlock::AppendProbabilities(const char *snp1, const char *snp2, Probabilities &prob1) {
	char key[1024];
	sprintf(key, "%sx%s", snp1, snp2);
	probLookup[key] = prob1;
}

void HaplotypeBlock::Append(size_t idx, Locus *snp) {
	if (idx > lastLocusIndex) {
		lastLocusIndex = idx;
		lastSnpLabel = snp->GetLabel();
	}
	if (idx < firstLocusIndex) {
		firstLocusIndex = idx;
		firstSnpLabel  = snp->GetLabel();
	}

	loci.push_back(HaplotypeBlock::LocusWrapper(idx, snp));
	double maf = snp->GetMinAlleleFreq();
	if (minAlleleFreq == NULL || minAlleleFreq->GetMinAlleleFreq() > maf)
		minAlleleFreq = snp;
	if (maxAlleleFreq == NULL || maxAlleleFreq->GetMinAlleleFreq() < maf)
		maxAlleleFreq = snp;
	totalMAF+=maf;
	assert(idx + 1 > lastIndex);
	lastIndex = idx;
}

bool HaplotypeBlock::operator<(const HaplotypeBlock& other) const {
	return loci.size() >= other.loci.size();
}

double HaplotypeBlock::GetAverageMAF() const {
	if (totalMAF > 0.0)
		return totalMAF/(double)loci.size();
	else 
		return 0.0;
}

size_t HaplotypeBlock::GetBlockDensity() const {
	return (GetLastLocation() - GetFirstLocation()) / BlockCount();
}

size_t HaplotypeBlock::GetFirstIndex() const {
	if (loci.size() > 0)	{
	
		vector<LocusWrapper>loci = this->loci;
	
		sort(loci.begin(), loci.end(), LTPos());
		return loci[0].index;
	}
	else	
		return 0;
}

size_t HaplotypeBlock::GetLastIndex() const {

	if (loci.size() > 0)	{
		vector<LocusWrapper>loci = this->loci;
	
		sort(loci.begin(), loci.end(), LTPos());
		return loci[loci.size() - 1].index;
	}
	else
		return 0;
}

uint HaplotypeBlock::GetFirstLocation() const {
	//OK, for ease of reading, we'll sort them- but we want to make a copy first
	vector<LocusWrapper>loci = this->loci;

	sort(loci.begin(), loci.end(), LTPos());

	return loci[0].locus->GetLocation();
}

uint HaplotypeBlock::GetLastLocation() const {
	//OK, for ease of reading, we'll sort them- but we want to make a copy first
	int lociCount = this->loci.size();
	vector<LocusWrapper>loci = this->loci;

	sort(loci.begin(), loci.end(), LTPos());

	return loci[lociCount-1].locus->GetLocation();
}

void HaplotypeBlock::Report(ostream& os) {
	size_t end=loci.size() - 1;
	os<<"Block: "<<loci[0].locus->GetLabel()<<" - "<<loci[end].locus->GetLabel()<<"  Size: "<<end+1<<" Average MAF: "<<GetAverageMAF()<<" BD: "<<GetBlockDensity()<<"\n";
}


static string tableHeaderColor = "#000000";
void HaplotypeBlock::HtmlSummaryBegin(ostream &os) {
	os<<"\n<DIV id=\"HaplotypeSummaryTable\"><!Summary Table for remaining blocks>\n";
	//os<<"<TABLE cellspacing=\"1\" cellpadding=\"1\" border=\"1\"><CAPTION>Block Summary</CAPTION>\n";
	os<<"<TABLE><CAPTION>Block Summary</CAPTION>\n";
	os<<"\t<tr><TH>First Snp</TH><TH>Last Snp</th><th>Block Size</th><th>Avg. Min. Allele Freq.</TH><TH>Block Density</TH></STRONG></tr>\n";
}
void HaplotypeBlock::HtmlSummary(ostream& os, const char *linkToDetails) {
	size_t end=loci.size() - 1;

	static string rowBackgrounds[] = {"\t<TR><TD>", "\t<TR class=\"invert\"><TD>"};

	int idx = rowsWritten++ %2;
	os<<rowBackgrounds[idx];
	if (linkToDetails) 
		os<<"<A HREF=#"<<linkToDetails<<">"<<loci[0].locus->GetLabel()<<"</A>";
	else
		os<<loci[0].locus->GetLabel();
	os<<"</TD><TD>"<<loci[end].locus->GetLabel()
		<<"</TD><TD>"<<end+1
		<<"</TD><TD>"<<GetAverageMAF()
		<<"</TD><TD>"<<GetBlockDensity()
		<<"</TD><TR>\n";
}

void HaplotypeBlock::HtmlSummaryEnd(ostream &os) {
	os<<"</TABLE></DIV>\n\n";
}

string HaplotypeBlock::HtmlDetailed(ostream &os, uint imgWidth, const char *dPrimeFilename, const char *rSquaredFilename, const char *dPrimeDetails, const char *rSquaredDetails) {
	size_t lociCount = this->loci.size();
	vector<LocusWrapper>loci = this->loci;
	string anchor=firstSnpLabel+"-"+lastSnpLabel;
	//Let's get rid of the "/"s in the filename

	size_t width = maxImageWidth;
	if (imgWidth < maxImageWidth)
		width = imgWidth;

	size_t offset = GetFirstIndex();

	sort(loci.begin(), loci.end(), LTMAF());

	os<<"<DIV class='spacer'></DIV>\n";
	os<<"<DIV class='spacer'></DIV>\n";
	os<<"\n<DIV id=\"HaplotypeBlock"<<firstSnpLabel<<"-"<<lastSnpLabel<<"\" width="<<width<<">\n";
	os<<"<H2><CENTER><A name="<<anchor<<">Haplotype Block "<<firstSnpLabel<<" to "<<lastSnpLabel<<"</A></CENTER></H2>\n";

	if (dPrimeFilename) {
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
		

		os<<"\t<DIV class='spacer'></DIV><DIV class='image-frame'><IMG SRC='"<<picFirst<<"' border=\"0\" alt=\"Haplotype Plot\" id=\""<<dPrimeFilename<<"\" width="<<width<<" onClick=\"window.open('"<<zoomedOut<<"')\"> </IMG></DIV>\n";
	}

	if (rSquaredFilename) {
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

		os<<"\t<DIV class='spacer'></DIV><DIV class='image-frame'><IMG SRC='"<<picFirst<<"' border=\"0\" alt=\"Haplotype Plot\" id=\""<<rSquaredFilename<<"\" width="<<width<<" onClick=\"window.open('"<<zoomedOut<<"')\"> </IMG></DIV>\n";
	}		
	

	os<<"\t<CENTER><TABLE><CAPTION>Block Details</CAPTION><TR><TH>Size</TH><TH>Block Density</TH><TH>Avg. MAF</TH></TR>\n";
	os<<"\t<TR class='invert'><TD>"<<lociCount<<"</TD> <TD>"<<GetBlockDensity()<<"</TD> <TD>"<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<GetAverageMAF()<<"</TD></TR>\n\t</TABLE>";

	os<<"\t<DIV class='spacer'></DIV><TABLE class=\"SummaryTable\"><CAPTION>Block Constituents</CAPTION>\n";
	os<<"\t\t<TR><TH>Index</TH><TH>Label</th><th>Location</th><th>Avg. Min. Allele Freq.</TH></tr>\n";

	for (size_t i=0; i<lociCount; i++) {
		LocusWrapper &lw = loci[i];
		Locus *loc = lw.locus;

		int idx = i%2;
		if (idx == 0)
			os<<"<TR class=\"invert\">";
		else
			os<<"<TR>";
	
    	os<<"<TD>"<<setw(12)<<right<<lw.index - offset + 1<<"("<<lw.index<<")"<<
			"</TD><TD>"<<setw(12)<<right<<loc->GetLabel()<<
			"</TD><TD>"<<setw(10)<<loc->GetLocation()<<
			"</TD><TD>"<<setw(12)<<setprecision(5)<<loc->GetMinAlleleFreq()<<"</TD></TR>\n";	
	}
	os<<"\t</TABLE>\n";
	os<<"<P><DIV class='spacer'></DIV><TABLE><CAPTION>Probability Components</CAPTION>"
	  <<"<TR><TH>Snps</TH><TH>cAB(pAB)</TH>"
	  <<"<TH>CAb(pAb)</TH><TH>caB(paB)</TH><TH>cab(pab)</TH><TH>maf(A)</TH><TH>maf(B)</TH></TR>\n";
	
	map<string, Probabilities>::iterator cur = probLookup.begin();
	map<string, Probabilities>::iterator end = probLookup.end();

	size_t i=0;	

	do {
		Probabilities &p = cur->second;
		int idx = ++i%2;
		if (idx == 0) {
			os<<"<TR class=\"invert\">";
		}
		else
			os<<"<TR>";
	
		os<<"<TD>"<<cur->first<<"</TD>"
		  <<"<TD>"<<p.cAB<<"("<<p.pAB()<<")</TD><TD>"<<p.cAb<<"("<<p.pAb()<<")</TD>"
		  <<"<TD>"<<p.caB<<"("<<p.paB()<<")</TD><TD>"<<p.cab<<"("<<p.pab()<<")</TD>"
		  <<"<TD>"<<p.mafA<<"</TD><TD>"<<p.mafB<<"</TD></TR>\n";
	}
	while (++cur != end);
	os<<"\t</TABLE>\n</DIV>\n";
	return anchor;
}


void HaplotypeBlock::DetailedReport(ostream &os) {
	//OK, for ease of reading, we'll sort them- but we want to make a copy first
	size_t lociCount = this->loci.size();
	vector<LocusWrapper>loci = this->loci;

	size_t offset = GetFirstIndex();

	sort(loci.begin(), loci.end(), LTMAF());

	os<<"\tHaplotype Block, Size: "<<lociCount<<" - "<<GetBlockDensity()<<"\n";
	os<<"\t\tBlock Constituents: "
			<<minAlleleFreq->GetLabel()<<" - "
			<<maxAlleleFreq->GetLabel()<<" Avg MAF: "
			<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)<<GetAverageMAF()<<"\n";
    os<<setw(12)<<right<<"Index"<<setw(12)<<right<<"Label"<<setw(10)<<"Location"<<setw(12)<<setprecision(5)<<"Min. Al. Freq"<<"\n";
	for (size_t i=0; i<lociCount; i++) {
		LocusWrapper &lw = loci[i];
		Locus *loc = lw.locus;
    	os<<setw(12)<<right<<lw.index - offset + 1<<
			setw(12)<<right<<loc->GetLabel()<<
			setw(10)<<loc->GetLocation()<<
			setw(12)<<setprecision(5)<<loc->GetMinAlleleFreq()<<"\n";	
	}
	

}

size_t HaplotypeBlock::BlockCount() const {
	return loci.size();
}

/**
 * @brief To calculate haplotype frequencies, we need to iterate over the contents of the pool and 
	      increment counts just like we did for ChromPool::CalculateFrequencies()
 */

}
