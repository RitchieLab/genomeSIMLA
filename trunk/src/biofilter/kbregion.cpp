//
// C++ Implementation: kbregion.cpp
//
// Description: Region based on knowledge based sources
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "kbregion.h"
#include <iomanip>


using namespace std;
using namespace Utility;

namespace Biofilter {
namespace Knowledge {
KbRegion::KbRegion() : KbEntity(0, "", ""),Region("", 0),
	startPosition(0), endPosition(0), origStart(0), origEnd(0), chromosome(""), snpLookup(NULL) { }

KbRegion::KbRegion(uint geneID, uint start, uint end, const char *chrom, const char *ensembl, const char *desc, SnpManager* snpLookup)
	: KbEntity(geneID, ensembl, desc), Region(ensembl, geneID), startPosition(start), endPosition(end),
		origStart(start), origEnd(end), chromosome(chrom), snpLookup(snpLookup) { }

KbRegion::KbRegion(const KbRegion& other):
	KbEntity(other), Region(other), startPosition(other.startPosition), endPosition(other.endPosition),
	origStart(other.origStart), origEnd(other.origEnd), chromosome(other.chromosome), snpLookup(other.snpLookup) { }

KbRegion::~KbRegion() { cerr<<"!!"<<Name()<<" is going away!\n"; }


bool KbRegion::AddSnp(uint snp) {
	bool success=false;

	if (snp >= startPosition && snp <= endPosition) {
		snps.insert(snp);
		success=true;
	}
	return success;
}


uint KbRegion::AssociateSNPs() {
	return snpLookup->GetSNPs(chromosome.c_str(), startPosition, endPosition, snps);
}

string KbRegion::RegionName() {
	return CommonName();
}
uint KbRegion::SnpCount(set<uint>& snpList) {
	set<uint> common;
	set_union(snps.begin(), snps.end(), snpList.begin(), snpList.end(), inserter(common, common.begin()));
	return common.size();
}

uint KbRegion::SnpCount() {
	return snps.size();
}

uint KbRegion::Start() {
	return startPosition;
}

uint KbRegion::OriginalStart() {
	return origStart;
}

uint KbRegion::OriginalEnd() {
	return origEnd;
}

uint KbRegion::End() {
	return endPosition;
}

KbRegion::SnpType::Type KbRegion::GetSnpTypeByRS(uint rsID) {
	SNPSet snpHits;
	if (snpManager->GetSNPs(rsID, snpHits) > 0) {
		if (snpHits.size() > 1) 
			return SnpType::Unknown;
		SNP_Details snp = snpManager->GetDetails(*snpHits.begin());
		uint pos = snp.position;
		if (pos > origStart && pos < origEnd)
			return SnpType::Interior;
		if (pos > startPosition && pos < endPosition)
			return SnpType::Flanking;
	}
	return SnpType::Exterior;
}

void KbRegion::GetBounds(uint& start, uint& end) {
	start = startPosition;
	end = endPosition;
}

void KbRegion::GetOrigBounds(uint& start, uint& end){
	start = origStart;
	end	= origEnd;
}

void KbRegion::SetBounds(uint start, uint end) {
	startPosition	= start;
	endPosition		= end;
}

void KbRegion::SetOrigBounds(uint start, uint end) {
	origStart		= start;
	origEnd			= end;
}

string KbRegion::Chromosome() {
	return chromosome;
}


int KbRegion::GetSnpCoverage(set<uint>& snpList, set<SNP_Details>& snpDetails) {
	uint count = 0;

	set<uint>::iterator itr = snps.begin();
	set<uint>::iterator end = snps.end();
	set<uint>::iterator notFound = snpList.end();
	set<uint> locals;
	int chromosome = ChromToInt(this->chromosome.c_str());
	while (itr != end) {
		SNP_Details snp = snpManager->GetDetails(*itr);
		assert(snp.chromosome==chromosome);
		if (snpList.find(*itr) != notFound) {
			count++;
			snpDetails.insert(snp);
		}
		itr++;
	}
	return count;
}

uint KbRegion::CollectSnpDetails(set<SNP_Details>& collection) {
	set<uint>::iterator itr = snps.begin();
	set<uint>::iterator end = snps.end();

	uint count = 0;
	while (itr != end) {
		uint pos = *itr++;
		collection.insert(snpLookup->GetDetails(pos));
		count++;
	}
	return count;
}

std::set<uint> KbRegion::SNPs() {
	std::set<uint> rsIDs;

	std::set<uint>::iterator itr = snps.begin();
	std::set<uint>::iterator end = snps.end();

	while (itr != end) {
		uint pos = *itr++;
		SNP_Details snp = snpLookup->GetDetails(pos);
		rsIDs.insert(snp.rsID);
	}
	return rsIDs;
}

set<uint> KbRegion::GetSnpCoverage(set<uint>& snpList) {
	set<uint>::iterator itr = snps.begin();
	set<uint>::iterator end = snps.end();
	set<uint>::iterator notFound = snpList.end();
	set<uint> locals;
	int chromosome = ChromToInt(this->chromosome.c_str());
	while (itr != end) {
		SNP_Details snp = snpManager->GetDetails(*itr);
		assert(snp.chromosome==chromosome);
		if (snpList.find(*itr) != notFound)
			locals.insert(snp.rsID);
		itr++;
	}
	return locals;
}

uint KbRegion::ListAssociations(int tabCount, ostream& os) {
	if (snps.size() > 0) {
		PrintTabs(tabCount, os);
		set<uint>::iterator itr = snps.begin();
		set<uint>::iterator end = snps.end();
		os<<CommonName()<<" ( ";
		while (itr != end) {
			SNP_Details snp = snpManager->GetDetails(*itr);
			os<<"rs"<<snp.rsID<<" ";
			itr++;
		}
		os<<" )\n";
	}
	return snps.size();
}

}
}
