//
// C++ Implementation: ldcorrection
//
// Description: 
//
//
// Author: Eric Torstenson <torstees@torstensonx.mc.vanderbilt.edu>, (C) Marylyn Ritchie 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ldcorrection.h"
#include <fstream>
#include <iomanip>
namespace Biofilter {

using namespace soci;
using namespace std;
using namespace Utility;

LdCorrection::LdCorrection()  { }

LdCorrection::~LdCorrection() { 
	PurgeRegions();
}


void LdCorrection::LoadConfiguration(soci::session& sociDB, const char *cfg) {
	ifstream ifile(cfg);
	if (!ifile.good()){ 
		cerr<<"Unable to read from file, cfg. Unable to load ld data.";
		exit(1);
	}
	int lineCount = 0;
	while (!ifile.eof()) {
		char line[4096];
		ifile.getline(line, 4096);
		if (line[0] != '#') {
			stringstream ss(line);
			if (lineCount++ == 0) 
				ss>>ldName>>ldComment;
			else {
				string cmd="";
				ss>>cmd;
				if (cmd == "RSQUARED")
					LoadValuesRS(sociDB, ss);
				else if (cmd == "DPRIME")
					LoadValuesDP(sociDB, ss);
				else
					if (cmd.length() > 0) 
						ldFilenames.push_back(line);
			}
		}
	}
}

void LdCorrection::LoadValuesRS(soci::session& sociDB, istream& input) {
	while (!input.eof()) {
		float val=0.0;
		input>>val;
		string pop = GetPopulationName("RS", val);
		string desc = GetPopDescription("RSquared", val);
		if (val > 0.0) {
			uint popID = InitPopulation(sociDB, pop.c_str(), desc.c_str());
			RegionSpline::AddRS(val, popID);
		}
	}
}

void LdCorrection::LoadValuesDP(soci::session& sociDB, istream& input) {
	while (!input.eof()){ 
		float val=0.0;
		input >> val;
		string pop = GetPopulationName("DP", val);
		string desc = GetPopDescription("DPrime", val);
		if (val > 0.0) {
			uint popID = InitPopulation(sociDB, pop.c_str(), desc.c_str());
			RegionSpline::AddDP(val, popID);
		}
	}
}

string LdCorrection::GetPopDescription(const char *stat, float threshold) {
	stringstream desc;
	desc<<ldName<<" Population. "<<stat<<" cutoff of "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(2)<<threshold;
	return desc.str();
}

string LdCorrection::GetPopulationName(const char *type, float threshold) {
	stringstream pop;
	pop<<ldName<<"-"<<type<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(2)<<threshold;
	//sociDB<<"SELECT population_id FROM populations WHERE population_label=:type", into(popID), use(string(type));
	return pop.str();
}

int LdCorrection::GetPopID(soci::session& sociDB, const char* type, float threshold) {
	string popName = GetPopulationName(type, threshold);
	sociDB<<"SELECT population_id FROM populations WHERE population_label=:type", into(popID), use(popName);
	return popID;
}


void LdCorrection::Process(soci::session& sociDB, const char *variationFilename) {
	vector<string>::iterator itr = ldFilenames.begin();
	vector<string>::iterator end = ldFilenames.end();

	while (itr != end)  {
		stringstream ss(*itr++);
		string chromosome="", filename = "";
		ss>>chromosome>>filename;
		int chrom = ChromToInt(chromosome.c_str());
		if (chrom>0) 
			ProcessLD(sociDB, chrom, filename.c_str(), variationFilename);
		else {
			cerr<<"Invalid LD configuration: LD Files must be preceded with the chromosome from which they are drawn from (1..22, X, Y, MT)\n";
			exit(1);
		}
	}
}

int LdCorrection::InitPopulation(soci::session& sociDB, const char *pop, const char *popDesc) {
	popID = -1;
	sociDB<<"SELECT population_id FROM populations WHERE population_label=:name", into(popID), use(string(pop));
	if (popID > 0) {
		cerr<<"Clearing out all bounds associated with population "<<popID<<"\n";
		sociDB<<"DELETE FROM region_bounds WHERE population_id=:id", use(popID);
		sociDB<<"DELETE FROM populations WHERE population_id=:id", use(popID);
	} else {
		sociDB<<"SELECT MAX(population_id) FROM populations", into(popID);
		popID++;
	}
	sociDB<<"INSERT INTO populations VALUES (:id, :name, :ldComment, :desc)", use(popID), use(string(pop)), use(string(popDesc)), use(string(ldComment));
	stringstream sql;
	sql <<"INSERT INTO region_bounds SELECT gene_id, "<<popID<<", start, end FROM region_bounds WHERE population_id=0";
	sociDB<<sql.str();
	return popID;
}


void LdCorrection::PurgeRegions() {
	RBTreeNode<uint, RegionSpline*> * node = regions.GetFirst();

	while (node) {
		delete node->GetData();
		node = node->GetNext();
	}
	regions.Clear();
}

void LdCorrection::LoadGenes(soci::session& sociDB, int chromosome) {
	rowset<row> genes = (sociDB.prepare<<"SELECT gene_id, chrom, start, end FROM regions NATURAL JOIN (SELECT gene_id, start, end FROM region_bounds WHERE population_id=0) b WHERE chrom = :chrom", use(IntToChrom(chromosome)));
	PurgeRegions();
	cerr<<"Collecting appropriate genes!\n";
	for (rowset<row>::const_iterator itr = genes.begin(); itr != genes.end(); ++itr) {
		uint chrom, geneID, start, end;		
		row const& row = *itr;
		geneID = row.get<int>(0);
		string chromosome = row.get<string>(1);
		chrom = ChromToInt(chromosome.c_str());
		start  = row.get<int>(2);
		end    = row.get<int>(3);
		regions.Add(start, new RegionSpline(geneID, chrom, start, end));
	}
	cerr<<regions.GetCount()<<"\n";

}


void LdCorrection::ProcessLD(soci::session& sociDB, int chrom, const char *ldSource, const char *variationFilename) {
	ifstream file(ldSource);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<ldSource<<"\n";
		exit(0);
	}

	//First, we need to collect all of the rsIDs that are found in the LD file
	set<uint> snps;
	uint lastID = 0;
	while (!file.eof() && file.good()) {
		string junk;
		string rs1, rs2;
		file>>junk>>junk>>junk>>rs1>>rs2>>junk>>junk>>junk>>junk;

		//We want to avoid reinserting the first one over and over, since we can very easily check that it's the same...#2 is more complicated, and the insert and search are about the same work
		rs1 = rs1.erase(0,2);
		rs2 = rs2.erase(0,2);
		uint curID = atoi(rs1.c_str());
		if (lastID != curID)
			snps.insert(curID);
		snps.insert(atoi(rs2.c_str()));
		lastID = curID;
	}
	file.clear();
	file.seekg(0, ios_base::beg);
	cerr<<snps.size()<<" SNPs found in file "<<ldSource<<".";cerr.flush();
	InitSNPs(chrom, snps, variationFilename);
	LoadGenes(sociDB, chrom);

	//If we have 0 regions, we can't really provide much help
	if (regions.GetCount() == 0)
		return;
	multimap<uint, uint>::iterator notFound = this->snps.end();

	while (!file.eof() && file.good()) {
		string pos1, pos2, pop, rs1, rs2;
		float dprime, rsquared, lod, fbin;
		file>>pos1>>pos2>>pop>>rs1>>rs2>>dprime>>rsquared>>lod>>fbin;
		if (pos1.length() == 0) 
			continue;
		
		rs1 = rs1.erase(0,2);
		rs2 = rs2.erase(0,2);

		uint s1 = atoi(rs1.c_str()), 
			s2 = atoi(rs2.c_str());

		//We are ignoring low LOD scores (unreliable signal) and missing SNPs (or those that have multiples)
		if (lod > 2.0 && (this->snps.count(s1) == 1 && this->snps.count(s2) == 1)) {
			uint left = this->snps.find(s1)->second;
			uint right = this->snps.find(s2)->second;
			RegionBoundariesNode *node = regions.GetFirst();		//regions.FindNearestMin(left);
			//This is to cover those cases where the gene extends past the first SNP (on the left) but the right side is still within the bounds...probably won't happen....
			if (node == NULL)
				node = regions.GetFirst();
			if (abs((int)left - (int)right) < 1000000) {
				while (node && node->GetKey() <= right) {
					RegionSpline *region = node->GetData();

					//To identify the SNP causing extreme boundary modifications

					if (region->AddSnps(left, right, chrom, dprime, rsquared))
						if ((left < region->start && abs((int)region->start - (int)left) > 600000 > 600000) || (right > region->end && abs((int)region->end - (int)right) > 600000)){
							cerr<<" -- rs"<<s1<<":"<<left<<" , rs:"<<s2<<":"<<right<<"   --  "<<region->start<<" , "<<region->end<<"\t\t"<<abs((int)region->start - (int)left)<<" - "<<abs((int)region->end - (int)right)<<"\n";
						}


					node = node->GetNext();
				}
			}
		}	
	}

	cerr<<"Committing the regions back to the database\n";
	RegionBoundariesNode *node = regions.GetFirst();
	//sociDB.begin();

	while (node) {
		node->GetData()->Commit(sociDB);
		node = node->GetNext();
	}
	//sociDB.commit();

/*
	std::stringstream ss;
	ss<<"BEGIN;\n";
	while (node) {
		node->GetData()->Commit(ss);
		if (ss.str().length() > 1000000) {
			ss<<"COMMIT;\n";
			sociDB<<ss.str();
			cout<<ss.str();
			ss.str("BEGIN;\n");
		}
		//node->GetData()->Commit(sociDB);
		node = node->GetNext();
	}
	if (ss.str().length() > 0) {
		ss<<"COMMIT;\n";

		sociDB<<ss.str();
		cout<<ss.str();
	}
	cout<<"Committing:";cout.flush();
	sociDB.commit();
	cout<<"...done!\n";
	 */
}



uint LdCorrection::InitSNPs(int chrom, std::set<uint>& snps, const char *fn) {
	Purge();
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

	cerr<<"**-> Snp Count("<<chrom<<") "<<this->snps.size()<<"\n";

	while (file.good()) {
		char label[3];
		file.read(label, 2);
		label[2]='\0';
		int chromosome = ChromToInt(label);
		int snpCount=0, maxPosition=0;
		file.read((char*)&snpCount, 4);
		file.read((char*)&maxPosition, 4);
		offset+= maxPosition;
		if (chrom != chromosome) {
			//Skip the rest of this chromosome
			file.seekg(snpCount*9, ios::cur);
		}
		else {
			Chromosome *newChrom = new Chromosome(label, offset);
			chrLookup[chromosome] = newChrom;
			posLookup[offset] = newChrom;
			if (file.good()) {
				cerr<<".";cerr.flush();
				for (int i=0; i<snpCount; i++) {
					int rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);
					if (chrom == chromosome) {
						if (snps.size() == 0 || (rs > 0 && snps.find(rs) != snps.end())) {
							idx=newChrom->AddSNP(pos, rs, role);
							this->snps.insert(pair<uint, uint>(rs, pos));
							count++;
						}
					}
				}
			}
			break;
		}
	}
	return count;
}


}
