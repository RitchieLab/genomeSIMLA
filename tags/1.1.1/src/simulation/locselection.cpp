//
// C++ Implementation: locusselection
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "locselection.h"

namespace Simulation {

namespace Visualization {

using namespace Utility;

float LocusSelection::mafWeight			= 100.0;
float LocusSelection::blockSizeWeight	= 50.0;
uint LocusSelection::maxReportedEntries	= 150;

void LocusSelection::SetDescription(const char *desc) {
	description = desc;
}

void LocusSelection::SetLabel(const char *lbl) {
	label = lbl;
}

bool LocusSelection::RemoveRange(size_t idx) {
	bool success = false;
	std::vector<SimpleRange<Locus> >::iterator itr = ranges.begin();
	size_t i=0;
	if (idx < ranges.size())  {
		while (i++ < idx)
			itr++;
		ranges.erase(itr);
		success=true;
	}
	return true;
}

string LocusSelection::GetConfiguration( const char *del) {
	stringstream ss;
	ss<<label<<del
	 	<<mafRange.target<<del
		<<mafRange.min<<del
		<<mafRange.max<<del
		<<blockSize.target<<del
		<<blockSize.min<<del
		<<blockSize.max<<del
		<<description;
	return ss.str();
}



		

bool LocusSelection::AddRange(Locus& start, Locus& stop) {
	ranges.push_back(SimpleRange<Locus>(start, stop));
	return true;
}

bool LocusSelection::SetRange(int idx, Locus &start, Locus &stop) {
	ranges[idx] = SimpleRange<Locus>(start, stop);
	return true;
}

bool LocusSelection::EvaluateRange(Locus &loc) {
	if (ranges.size() == 0) 
		return true;
	bool withinRange = false;

	vector<SimpleRange<Locus> >::iterator itr = ranges.begin();
	vector<SimpleRange<Locus> >::iterator end = ranges.end();

	while (!withinRange && itr != end) {
		withinRange = itr->Evaluate(loc) * mafWeight;
		itr++;
	}
	
	return withinRange;
}

void LocusSelection::ExtractRanges( LocusSelection &other) {
	std::vector<SimpleRange<Locus> >::iterator itr = other.ranges.begin();
	std::vector<SimpleRange<Locus> >::iterator end = other.ranges.end();

	while (itr != end) 
		ranges.push_back(*itr++);	
}

float LocusSelection::Evaluate(LocusReportEntry* locReport) {
	locReport->score = mafRange.Evaluate(locReport->locus.GetMinAlleleFreq()) * mafWeight;

	if (locReport->score > 0.0) {
		vector<BlockListNode *>::iterator itr = locReport->associatedBlocks.begin();
		vector<BlockListNode *>::iterator end = locReport->associatedBlocks.end();

		float blockScore = 0.0;
		while (itr != end) {
			uint snpCount = (*itr++)->GetBlockSize();
			float curBlockScore = blockSize.Evaluate(snpCount) * blockSizeWeight;
			if (curBlockScore > blockScore)
				blockScore = curBlockScore;
		}
		if (blockScore > 0.0) {
			if (blockSize.max > 0.0)
				locReport->score += blockScore;
		}
		else
			if (blockSize.min == 0 && locReport->associatedBlocks.size() == 0)
				locReport->score+=blockSizeWeight;
			else
				locReport->score = 0.0;

	}

	
	if (locReport->score > 0.0) {
		contents.Insert(locReport->score, locReport);
	}
	

	return locReport->score;
}

void LocusReportEntry::WriteToText(ostream &file) {
	file<<locus.GetLabel()<<" "<<score<<" "<<locus.GetLabel()<<" "<<locus.GetMinAlleleFreq()<<" "<<chromosomeID<<" ";
	for (uint n=0; n<associatedBlocks.size(); n++) {
		BlockListNode *block = associatedBlocks[n];
		file<<blockDescriptionFile<<" ("<<block->GetLabel()<<", "<<block->GetBlockSize()<<") ";
	}
	file<<"\n";
}

void LocusReportEntry::WriteHtmlHeader(ostream& report) {
	report<<"\t<TR><TH>Locus</TH><TH>Score</TH><TH>Freq (Allele 1)</TH><TH>Freq (Allele 2)</TH><TH>Chromosome</TH><TH>Associated Blocks</TH></TR>\n";
}

void LocusReportEntry::WriteHtmlReport(ostream& report, bool useAltColor) {
	float freqA = locus.Freq1();
	float freqa = locus.Freq2();

	if (freqa > freqA) 
		report<<"\t<TR class='invert'>";
	else
		report<<"\t<TR               >";
	report<<"<TD>"<<locus.GetLabel()<<" ("<<locus.GetChromID()+1<<":"<<locus.GetID()+1<<")</TD><TD>"
		<<score<<"</TD>";
	report<<"<TD>"<<freqA<<"</TD>";
	report<<"<TD>"<<freqa<<"</TD>";
	report<<"<TD>"<<chromosomeID<<"</TD>";
	report<<"<TD><TABLE><TR>";
	for (uint n=0; n<associatedBlocks.size(); n++) {
		BlockListNode *block = associatedBlocks[n];
		report<<"<TR><TD><A HREF='"<<blockDescriptionFile<<"#"<<block->GetLabel()<<"'>"<<block->GetLabel()<<" ("<<block->GetBlockSize()<<")</A></TD></TR>\n";
	}
	report<<"\t</TR></TABLE>\n";

}

void LocusSelection::WriteToText(ostream& report) {
	RBTreeNode<float, LocusReportEntry*> *cur = contents.GetLast();

	report<<"Locus Selection: "<<label<<"\n"<<description<<"\n";
	report<<"Target Freq, Min Freq, Max Freq, Target Bl Size, Min, Max Bl Size, Region Constraints\n";
	report<<mafRange.target<<", "<<mafRange.min<<", "<<mafRange.max<<",";
	report<<blockSize.target<<", "<<blockSize.min<<", "<<blockSize.max<<", ";

	std::vector<SimpleRange<Locus> >::iterator itr = ranges.begin();
	std::vector<SimpleRange<Locus> >::iterator end = ranges.end();

	if (itr == end) {
		report<<"no constraints";
	}
	cout<<"\n";
	while (itr != end) {
		report<<"\t"<<itr->min.GetLabel()<<" - "<<itr->max.GetLabel();
		itr++;
	}
	report<<"\n";

	if (cur) {	
		report<<"Matching Loci\n";
	
		for (uint i=0; cur && i<maxReportedEntries; i++) {
			LocusReportEntry *entry = cur->GetData();

			if (i==0)
				entry->WriteToText(report);

			assert(cur != cur->GetPrev());
			cur = cur->GetPrev();
		}
	
		report<<"\n";
	}
	else {
		report<<"---- none found\n";
	}
}

void LocusSelection::WriteHtmlReport(ostream& report) {
	RBTreeNode<float, LocusReportEntry*> *cur = contents.GetLast();
	report<<"<H2>Locus Selection: "<<label<<"</H2>\n"<<description<<"\n<P>";
	report<<"\n\t<P><TABLE><CAPTION>Locus Search Results</CAPTION>";
	report<<"<TR><TH colspan='3'>Minor Allele<P>Frequency</TH><TH colspan='3'>Block Size</TH><TH>Region </TH></TR>\n";
	report<<"<TR><TH>Target</TH><TH>Min</TH><TH>Max</TH><TH>Target</TH><TH>Min</TH><TH>Max</TH><TH>Constraints</TH></TR>\n";
	report<<"<TR name='Locus Search Results'><TD>"<<mafRange.target<<"</TD><TD>"<<mafRange.min<<"</TD><TD>"<<mafRange.max<<"</TD>";
	report<<"<TD>"<<blockSize.target<<"</TD><TD>"<<blockSize.min<<"</TD><TD>"<<blockSize.max<<"</TD>\n";

	std::vector<SimpleRange<Locus> >::iterator itr = ranges.begin();
	std::vector<SimpleRange<Locus> >::iterator end = ranges.end();
	report<<"<TD><TABLE>";
	if (itr == end) {
		report<<"<TR><TD>no constraints</TD></TR>\n";
	}
	while (itr != end) {
		report<<"<TR><TD>"<<itr->min.GetLabel()<<" - "<<itr->max.GetLabel()<<"</TD></TR>";
		itr++;
	}
	report<<"</TABLE></TD>\n</TABLE>\n<DIV class='spacer'></DIV>\n";

	if (cur) {	
		report<<"<TABLE name='"<<label<<"'><CAPTION>Matching Loci ("<<contents.GetCount()<<")</CAPTION>\n";
	
		bool alt=true;


		for (uint i=0; cur && i<maxReportedEntries; i++) {
			LocusReportEntry *entry = cur->GetData();

			if (i==0)
				entry->WriteHtmlHeader(report);

			alt = !alt;
			entry->WriteHtmlReport(report, alt);

			assert(cur != cur->GetPrev());
			cur = cur->GetPrev();
		}
	
		report<<"</TABLE>\n";
	}
	else {
		report<<"---- none found\n";
	}
	report<<"<DIV class='spacer'></DIV><DIV class='spacer'></DIV>\n";
}

LocusReportEntry::~LocusReportEntry() {
	BlockList::iterator itr = associatedBlocks.begin();
	BlockList::iterator end = associatedBlocks.end();

	while (itr != end) {
		delete (*itr++);
	}
}

void LocusSelection::Reset() {
	contents.Clear();
}


string LocusSelection::GetDescMAF() {
	return mafRange.GetDescription();
}


string LocusSelection::GetDescBlockSize() {
	return blockSize.GetDescription();
}

string LocusSelection::GetDescRanges() {
	stringstream ss;
	
	std::vector<SimpleRange<Locus> >::iterator itr = ranges.begin();
	std::vector<SimpleRange<Locus> >::iterator end = ranges.end();

	if (itr == end) {
		ss<<"no constraints";
	}
	while (itr != end) {
		ss<<itr->min.GetLabel()<<".."<<itr->max.GetLabel()<<" ";
		itr++;
	}

	return ss.str();
}

void LocusSelection::SetBlockSizeRange(LinearRange<uint> range) {
	blockSize = range;
}

void LocusSelection::SetMAFRange(LinearRange<float> range) {
	mafRange = range;
}



}

}
