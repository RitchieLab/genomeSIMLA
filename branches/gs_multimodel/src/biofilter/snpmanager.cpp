//
// C++ Implementation: snpmanager
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snpmanager.h"
#include "kbentity.h"
#include <fstream>
#include <iostream>
//#include <c++/4.0.0/bits/stl_bvector.h>

using namespace std;
using namespace Utility;

namespace Biofilter {

std::map<uint, std::string> Chromosome::roleDescription;

Chromosome::Chromosome(const char *label, uint offset) : offset(offset), label(label) { }
Chromosome::~Chromosome() { }
int Chromosome::CollectSnps(uint left, uint right, SNPSet& bag) {
	std::map<uint, uint>::iterator itr = snps.lower_bound(left);
	std::map<uint, uint>::iterator end = snps.upper_bound(right);
	int count = 0;
	while (itr != end) {
		count++; 
		bag.insert(itr++->first + offset);
	}

	return count;
}
SNP_Details Chromosome::GetDetails(uint pos) {
	SNP_Details details;
	if (snps.find(pos - offset) != snps.end()) {
			details.chromosome = Utility::ChromToInt(label.c_str());
			details.position = pos - offset;
			details.rsID = snps[details.position];
	}
	else
		cerr<<"\n?? "<<label<<" "<<pos<<" ( "<<pos-offset<<" )\n";
	return details;
}

void Chromosome::GetRsIDs(set<uint>& positions, set<uint>& rsIDs) {
	set<uint>::iterator itr = positions.begin();
	set<uint>::iterator end = positions.end();

	while (itr != end) {
		uint pos = *itr;
		if (snps.find(pos) != snps.end())
			rsIDs.insert(snps[pos]);
		itr++;
	}
}


bool Chromosome::operator<(const Chromosome& other) const {
	return offset < other.offset;
}
uint Chromosome::AddSNP(uint position, uint rsID, uint role) {
//cerr<<"+ "<<label<<" ( "<<offset<<" ) "<<position<<" rs"<<rsID<<" ( "<<position+offset<<" )\n";
	snps[position] = rsID;
	roles[position] = role;
	return position+offset;
}

uint Chromosome::GetRoleID(uint position) {
	if (roles.find(position) != roles.end())
		return roles[position];
	return 0;
}
void Chromosome::WriteMarkerInfo(ostream& os, bool detailed) {
	std::map<uint, uint>::iterator itr = snps.begin();
	std::map<uint, uint>::iterator end = snps.end();
	
	while (itr != end) {
		os<<"rs"<<itr->second<<"\t"<<itr->first<<"\t"<<label;
		if (detailed) {
			uint roleID = roles[itr->first];
			if (roleDescription.find(roleID) != roleDescription.end())
				os<<"\t"<<roleDescription[roleID];
			else
				os<<"\t"<<roleID;
		}
		os<<"\n";
		itr++;
	}
}

void Chromosome::PrintSNPs(ostream& os) {
	std::map<uint, uint>::iterator itr = snps.begin();
	std::map<uint, uint>::iterator end = snps.end();
	
	//This function doesn't work right yet. I don't have the gene information stored that way right now
	while (itr != end) {
		os<<"?@?@?@?@?@?@"<<"rs"<<itr->second<<"\t"<<itr->first<<"\t"<<label<<"\n";
		itr++;
	}
}





SnpManager::SnpManager() : filename("variations.bn") {	}
SnpManager::~SnpManager() { 
	Purge();
}


void SnpManager::Purge() {
	std::map<uint, Chromosome*>::iterator itr = posLookup.begin();
	std::map<uint, Chromosome*>::iterator end = posLookup.end();
	while (itr != end) {
		delete itr++->second;
	}
	posLookup.clear();
	chrLookup.clear();
	snps.clear();
}

int SnpManager::GetSNPs(const char *chromLabel, uint start, uint stop, SNPSet& snps) {
	int count = 0;
	int chrom = ChromToInt(chromLabel);
	if (chrLookup.find(chrom) != chrLookup.end()) 
		count = chrLookup[chrom]->CollectSnps(start, stop, snps);
	return count;
}
void SnpManager::GetDetails(SNPSet& snps, SnpDetailsCollection& details) {
	SNPSet::iterator itr = snps.begin();
	SNPSet::iterator end = snps.end();

	while (itr != end )
		details.insert(GetDetails(*itr++));
}
SNP_Details SnpManager::GetDetails(uint pos) {
	assert(posLookup.lower_bound(pos) != posLookup.end());
	return posLookup.lower_bound(pos)->second->GetDetails(pos);
}


uint SnpManager::InitSNPs(std::vector<uint>& snpSource, int chromosome, const char *filename, const char *mergedReport, bool doReportMerged) {
	uint convertedSNPs					= 0;
	stringstream ss;
	ifstream file(filename, ios::binary);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}
	set<uint>snps;
	snps.insert(snpSource.begin(), snpSource.end());
	uint count = 0;
	uint offset = 0;
	uint idx = 0;
	file.read((char*)&varVersion, 4);

	if (varVersion != varVersion) {
		cerr<<"Mismatched variation file. Expected: "<<varVersion<<" found "<<varVersion<<"\n";
	}
	while (file.good()) {
		char label[3];
		file.read(label, 2);
		label[2]='\0';
		int snpCount=0, maxPosition=0;
		file.read((char*)&snpCount, 4);
		file.read((char*)&maxPosition, 4);

		int chrom = ChromToInt(label);
		if (chrom != chromosome) {
			//Skip the rest of this chromosome
			file.seekg(snpCount*8, ios::cur);
		}
		else {
			Chromosome *newChrom = new Chromosome(label, offset);
			offset+= maxPosition;
			chrLookup[chrom] = newChrom;
			posLookup[offset] = newChrom;
			if (file.good()) {
				//cerr<<".";cerr.flush();
				for (int i=0; i<snpCount; i++) {
					int rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);
					if (rs > 0) {
						if (snps.size() == 0 || snps.find(rs) != snps.end()) {
							idx=newChrom->AddSNP(pos, rs, role);
							this->snps.insert(pair<uint, uint>(rs, idx));
							count++;
						} else {
							if (rsConversion.find(rs) != rsConversion.end()) {
								uint newrs					= rsConversion[rs];
								if (snps.find(newrs) != snps.end()) {
									idx = newChrom->AddSNP(pos, newrs, role);
									this->snps.insert(pair<uint, uint>(newrs, idx));
									ss<<label<<"\t"<<rs<<"\t"<<pos<<"\t"<<newrs<<"\n";
									convertedSNPs++;
									count++;
								}
							}
						}
					}
				}
			}
			//Break out of the while loop, since we've successfully loaded the contents of the chromsome
			break;
		}
	}
	Knowledge::KbEntity::snpManager = this;
	return count;
}


uint SnpManager::InitSNPs(std::vector<uint>& snpSource, const char *fn, const char *mergedReport, bool doReportMerged) {
	uint convertedSNPs					= 0;
	stringstream ss;
	set<uint> snps;

	snps.insert(snpSource.begin(), snpSource.end());

	if (fn)
		filename = fn;
	ifstream file(filename.c_str(), ios::binary);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}
	uint count = 0;
	uint offset = 0;
	uint idx = 0;
	file.read((char*)&varVersion, 4);

	if (varVersion != varVersion) {
		cerr<<"Mismatched variation file. Expected: "<<varVersion<<" found "<<varVersion<<"\n";
	}
	while (file.good()) {
		char label[3];
		file.read(label, 2);
		label[2]='\0';
		int snpCount=0, maxPosition=0;
		file.read((char*)&snpCount, 4);
		file.read((char*)&maxPosition, 4);

		int chrom = ChromToInt(label);
		Chromosome *newChrom = new Chromosome(label, offset);
		offset+= maxPosition;
		chrLookup[chrom] = newChrom;
		posLookup[offset] = newChrom;
		if (file.good()) {
			//cerr<<".";cerr.flush();
			for (int i=0; i<snpCount; i++) {
				int rs=0, pos=0, role=0;
				file.read((char*)&rs, 4);
				file.read((char*)&pos, 4);
				file.read((char*)&role, 1);
				cerr<<"-- "<<rs<<" "<<pos<<" "<<role<<"\n";
				if (rs > 0) {
					if (snps.size() == 0 || snps.find(rs) != snps.end()) {
						idx=newChrom->AddSNP(pos, rs, role);
						this->snps.insert(pair<uint, uint>(rs, idx));
						count++;
					} else {
						if (rsConversion.find(rs) != rsConversion.end()) {
							uint newrs					= rsConversion[rs];
							if (snps.find(newrs) != snps.end()) {
								idx = newChrom->AddSNP(pos, newrs, role);
								this->snps.insert(pair<uint, uint>(newrs, idx));
								ss<<label<<"\t"<<rs<<"\t"<<pos<<"\t"<<newrs<<"\n";
								convertedSNPs++;
								count++;
							}
						}
					}
				}
			}
		}
	}
	if (doReportMerged) {
		cerr<<"SNPs Converted (rs_merge): "<<convertedSNPs<<"\n";
		ofstream file(mergedReport);
		file<<convertedSNPs<<"Chrom\tNew RS\tPosition\tOld RS:\n"<<ss;
		file.close();
	}
	Knowledge::KbEntity::snpManager = this;
	assert(this->snps.size() == count);
	return count;
}

uint SnpManager::InitSNPs(std::set<uint>& snps, const char *fn, const char *mergedReport, bool doReportMerged) {
	uint  convertedSNPs							= 0;
	stringstream ss;
	if (fn)
		filename										= fn;
	ifstream file(filename.c_str(), ios::binary);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}
	uint count										= 0;
	uint offset										= 0;
	uint idx											= 0;

	set<uint> foundSnps;
	//uint varVersion = 0;
	file.read((char*)&varVersion, 4);

	if (varVersion != varVersion) {
		cerr<<"Mismatched variation file. Expected: "<<varVersion<<" found "<<varVersion<<"\n";
	}
	while (file.good()) {
		char label[3];
		file.read(label, 2);
		label[2]='\0';
		int snpCount								=0,
			 maxPosition							=0;
		file.read((char*)&snpCount, 4);
		file.read((char*)&maxPosition, 4);

		int chrom = ChromToInt(label);
		Chromosome *newChrom = new Chromosome(label, offset);
		offset+= maxPosition;
		chrLookup[chrom] = newChrom;
		posLookup[offset] = newChrom;
		if (file.good()) {
			//cerr<<".";cerr.flush();
			for (int i=0; i<snpCount; i++) {
				int rs=0, pos=0, role=0;
				file.read((char*)&rs, 4);
				file.read((char*)&pos, 4);
				file.read((char*)&role, 1);

				if (rs > 0) {
					if (snps.size() == 0 || snps.find(rs) != snps.end()) {
						idx=newChrom->AddSNP(pos, rs, role);
						this->snps.insert(pair<uint, uint>(rs, idx));
						foundSnps.insert(rs);
						count++;
					} else {
						if (rsConversion.find(rs) != rsConversion.end()) {
							uint newrs					= rsConversion[rs];
							if (snps.find(newrs) != snps.end()) {
								cerr<<"-- Converted SNP: "<<rs<<" -> "<<newrs<<"\t"<<rsConversion.size()<<"\n";
								idx = newChrom->AddSNP(pos, newrs, role);
								this->snps.insert(pair<uint, uint>(newrs, idx));
								ss<<label<<"\t"<<rs<<"\t"<<pos<<"\t"<<newrs<<"\n";
								convertedSNPs++;
								count++;
								foundSnps.insert(newrs);
							}
						}
					}
				}
			}
		}
	}


	if (doReportMerged) {
		set<uint> missingSNPs;
		
		set_difference(snps.begin(), snps.end(), foundSnps.begin(), foundSnps.end(), inserter(missingSNPs, missingSNPs.begin()));
	
		cerr<<"SNPs Converted (rs_merge): "<<convertedSNPs<<"\n";
		cerr<<"SNPs remaining unconverted: "<<missingSNPs.size()<<"\n";

		set<uint>::iterator itr = missingSNPs.begin();
		set<uint>::iterator end = missingSNPs.end();
		ofstream bs("bad-snps.csv");

		while (itr != end)  {
			bs<<"rs"<<*itr++<<"\n";
		}
 
		ofstream file(mergedReport);
		file<<convertedSNPs<<"Chrom\tNew RS\tPosition\tOld RS:\n"<<ss.str();
		file.close();
	}


	Knowledge::KbEntity::snpManager = this;
	assert(this->snps.size() == count);
	return count;
}

void SnpManager::WriteMarkerInfo(ostream& os, bool detailed) {
	std::map<uint, Chromosome*>::iterator itr = posLookup.begin();
	std::map<uint, Chromosome*>::iterator end = posLookup.end();

	while (itr != end) 
		itr++->second->WriteMarkerInfo(os, detailed);
}

void SnpManager::PrintSNPs(ostream& os) {
	std::map<uint, Chromosome*>::iterator itr = posLookup.begin();
	std::map<uint, Chromosome*>::iterator end = posLookup.end();

	while (itr != end) 
		itr++->second->PrintSNPs(os);
}



/**
 * It should be noted that this returns the position relative to the beginning of the
 * genome....not the snps position on the chromosome...
 */
int SnpManager::GetSNPs(uint rsID, SNPSet& stache) {
	int count = 0;
	std::multimap<uint, uint>::iterator itr = snps.lower_bound(rsID);
	if (itr != snps.end()) {
		std::multimap<uint, uint>::iterator end = snps.upper_bound(rsID);
		while (itr != end) {
			count++; 
			uint pos = (itr++)->second;
			stache.insert(pos);
		}
	}
	return count;
}

}
