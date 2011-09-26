#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <mysql.h>
#include <iostream>
#include <sstream>
//#include <c++/4.0.0/bits/stl_bvector.h>
using namespace std;
MYSQL *connection, mysql;
MYSQL_RES *result;
MYSQL_ROW row;


using namespace std;

struct Locus {
	string rsID;
	int position;

	Locus() : rsID(""), position(-1) { }
	Locus(const char* rsid, int position) : rsID(rsid), position(position) { }
	Locus(const Locus& other) : rsID(other.rsID), position(other.position) {}
};

class DistanceCalculator {
public:
	DistanceCalculator() {}
	~DistanceCalculator() {}


	void Load(const char *filename, int start, int stop);

	float GetDistance(int bp);

protected:
	map<int, float> distances;
};

float DistanceCalculator::GetDistance(int bp) {
	map<int, float>::iterator itr = distances.find(bp);

	//If we don't find the bp locaton in our map, we will extrapolate the 
	//distance using a weighted average of the two surrounding ones..
	if (itr == distances.end()) {
		//If we crash here, it's likely because we fell off of the tree...ouch!
		int pLoc=0, nLoc=0;
		float pDist=0.0, nDist=0.0;

		itr = distances.lower_bound(bp);
		nLoc = itr->first;
		nDist = itr->second;

		itr--;
		pLoc = itr->first;
		pDist = itr->second;

		float length = nLoc;
		if (pLoc > nLoc)
			length = nLoc-pLoc;
		
		
		float prop = float(bp-pLoc) / length;
		cout<<"--gd ( "<<bp<<" ) = "<<nLoc<<"\t"<<nDist<<"\t"<<length<<"\t"<<prop<<"\t"<<prop*pDist +(1.0-prop)*nDist<<"\n";
		return prop*pDist + (1.0-prop)*nDist;
	}
	return itr->second;
}

void DistanceCalculator::Load(const char *filename, int start, int stop) {
	ifstream file(filename);
	cerr<<"Loading Distances from file, "<<filename<<"\n";
	float lastPosition = 0.0;

	char header[1024];
	file.getline(header, 1024);
	while (file.good() && !file.eof()) {
		file.getline(header, 1024);
		stringstream ss(header);
		int pos = 0;
		float cRate = 0.0;
		float gmap = 0.0;

		ss>>pos>>cRate>>gmap;
/*		if (pos > start && pos < stop) {
			cout<<pos<<"\t"<<cRate<<"\t"<<gmap<<"\n";
		} */
		if (pos>0) {
			distances[pos] = gmap - lastPosition;
			lastPosition = gmap;
		}

	}
}

float GetFreq() {
	return 0.5;
}

class Individual {
public:
	Individual(ostream* file, map<int, Locus>& snps) : id(""), os(file), completedMap(false), snps(snps) { }
	~Individual() { }

	void Reset(const char* id);
	bool AddGenotype(const char* id, int position, const char *genotype);
	void WriteToPed();
	void WriteMapFile(ostream& file);
protected:
	string id;
	ostream *os;
	vector<int> positions;
	map<int, string> genotypes;
	map<int, string> genotypes2;
	bool completedMap;
	map<int, Locus> snps;
};

void Individual::Reset(const char* id) {
	if (genotypes.size() > 0) {
		genotypes.clear();
		genotypes2.clear();
		completedMap = true;
	}
	this->id = id;
}

void Individual::WriteMapFile(ostream& file) {
	vector<int>::iterator itr = positions.begin();
	vector<int>::iterator end = positions.end();
	cerr<<"Total Locus Count: "<<positions.size()<<"\n";
	while (itr != end) {
		file<<snps[*itr].rsID<<" "<<snps[*itr].position<<"\n";
		itr++;
	}
}

bool Individual::AddGenotype(const char* id, int position, const char* genotype) {

	if (this->id != string(id)) {
		WriteToPed();
		Reset(id);
	}

	if (genotypes.find(position) != genotypes.end())
		genotypes[position] = genotype;
	else
		genotypes2[position] = genotype;
	if (!completedMap)
		positions.push_back(position);
	return false;
}

void Individual::WriteToPed() {
	if (genotypes.size() > 0) {
		vector<int>::iterator itr = positions.begin();
		vector<int>::iterator end = positions.end();
		*os<<id<<" "<<id<<" 0 0 0 0";
		while (itr != end) {
			*os<<" "<<genotypes[*itr]<<" "<<genotypes2[*itr];
			itr++;
		}
		*os<<"\n";
	}
}

void ShowEthnicities(MYSQL* cnx) {
	char sql[1024];
	sprintf(sql, "SELECT race FROM hapmap.1000genomes_race");
	mysql_query(cnx, sql);
	MYSQL_RES *result = mysql_store_result(cnx);
	while ((row=mysql_fetch_row(result))!=NULL) {
		cerr<<"\t"<<row[0]<<"\n";
	}

}

int GetEthMap(MYSQL* cnx, const char* eth) {
	char sql[1024];
	sprintf(sql, "SELECT id FROM hapmap.1000genomes_race WHERE race='%s'", eth);
	//SELECT id, race FROM hapmap.1000genomes_race");
	int query_state = mysql_query(cnx, sql);
	if (query_state != 0) {
		cerr<<mysql_error(cnx)<<"\n";
		abort();
	}
	MYSQL_RES *result = mysql_store_result(cnx);
	if ((row=mysql_fetch_row(result))!=NULL)
		return atoi(row[0]);

	cerr<<"Unrecognized ethnicity: "<<eth<<"\n";
	ShowEthnicities(cnx);
	return -1;
}

int main(int argc, char**argv) {
	//Connect to mysql-right now, we are hard coding it to justin's account. Woohoo!
	mysql_init(&mysql);
	connection = mysql_real_connect(&mysql, "munster", "justin", "lambchop@227","hapmap", 0, 0, 0);
	if (argc == 6) {
		string eth = argv[1];
		string chr = argv[2];
		int begin = atoi(argv[3]);
		int end   = atoi(argv[4]);
		DistanceCalculator dc;
		map<int, Locus> snps;

		
		uint ethID = GetEthMap(connection, eth.c_str());
		if (connection) {
			char sql[1024];
			sprintf(sql, "SELECT a.chr, a.var_id, a.pos FROM 1000genomestest_map a WHERE a.chr=%s AND race_id=%d AND a.pos BETWEEN %d AND %d ORDER BY a.pos", chr.c_str(), ethID, begin, end);
			int query_state = mysql_query(connection, sql);
			if (query_state != 0) {
				cerr<<mysql_error(connection)<<"\n";
				return 1;
			}
			
			//Now that we have gotten a query set up, let's go ahead and load the genetic distance data from hapmap
			dc.Load(argv[5], begin, end);

			string nameBase = "chrom_" + chr + "_" + argv[3] + "_" + argv[4];
			ofstream locFile((nameBase + ".locus").c_str());

			locFile<<"Chr "<<chr<<" based on 1000 genomes with genetic distances calculated from hapmap distance data\n";
			locFile<<"\nLabel\tFreq Al1\tFreq Al2\tMap Dist.\tPosition\tDescription\n";
			result = mysql_store_result(connection);
			while ((row=mysql_fetch_row(result))!= NULL) {
				float freq = GetFreq();
				locFile<<row[1]<<"\t"<<freq<<"\t"<<1.0-freq<<"\t"<<dc.GetDistance(atoi(row[2]))<<"\t"<<row[2]<<"\tchr"<<row[0]<<"\n";
				snps[atoi(row[2])] = Locus(row[1], atoi(row[2]));
			}
			locFile.close();
			cout<<"GenomeSIMLA compatible locus file: "<<nameBase<<".locus\n";
			sprintf(sql, "SELECT a.chr, a.pos, a.map_loc FROM 1000genomestest_map a WHERE chr=%s AND race_id=%d AND pos BETWEEN %d AND %d ORDER BY a.pos", chr.c_str(), ethID, begin, end);
			cerr<<sql<<"\n";
			//sprintf(sql, "SELECT * FROM (SELECT a.chr, a.var_id, a.pos, b.indiv_id, group_concat(substr(b.genome, a.map_loc, 1) SEPARATOR ' ') FROM 1000genomestest_map a inner join 1000genomestest b on a.chr = b.chr WHERE b.race=%d AND a.chr=%s AND a.pos BETWEEN %d AND %d GROUP BY b.indiv_id, a.pos) m ORDER BY pos", ethID, chr.c_str(), begin, end);
	
			mysql_query(connection, sql);
			result = mysql_store_result(connection);
			map<int, int> positions;
			int first = -1, last = -1;
			while ((row=mysql_fetch_row(result)) != NULL) {
				last = atoi(row[2]);
				if (first < 0)
					 first = last;
				positions[last] = atoi(row[1]) - 1;
			}
			//sprintf(sql, "SELECT chr, indiv_id, substr(genome, %d, %d) FROM 1000genomestest WHERE race=%d AND chr=%s", first, last-first, ethID, chr.c_str());
			sprintf(sql, "SELECT chr, indiv_id, genome FROM 1000genomestest WHERE race=%d AND chr=%s", ethID, chr.c_str());
			cerr<<"\nTotal Positions: "<<positions.size()<<"\n"<<sql<<"\n";
			ofstream pedfile((nameBase + ".ped").c_str());
			Individual ind(&pedfile, snps);
			if (mysql_query(connection, sql) == 0) {
				result = mysql_store_result(connection);
				while ((row=mysql_fetch_row(result)) != NULL) {
					string genotypes = row[2];
					int genoCount = genotypes.size();
					string indID = row[1];

					cerr<<" -- "<<row[0]<<"\t"<<row[1]<<" -- \t"<<indID<<"\t"<<"\t"<<genotypes.size()<<"\n";
					
					map<int, int>::iterator pitr = positions.begin();
					map<int, int>::iterator pend = positions.end();

					while (pitr!=pend) {
						int idx = pitr->first;
						if (idx < genoCount) {
							char genotype = genotypes[idx];
							ind.AddGenotype(indID.c_str(), idx, &genotype);
						}
						pitr++;
					}
				}
				ind.WriteToPed();

			} else 
				cerr<<"Unable to perform query: "<<sql<<"\n";
			cout<<"Pedigree Data: "<<nameBase<<".ped\n";
			pedfile.close();
			ofstream mapFile((nameBase + ".map").c_str());
			ind.WriteMapFile(mapFile);
			mapFile.close();
			cout<<"Map File: "<<nameBase<<".map\n";
		}
		else {
			cerr<<mysql_error(&mysql)<<"\n";
			cerr<<"Unable to connect to the munster database for some reason.....grrr!\n";
		}


	}
	else {
		cerr<<"Usage: eth chrom bp-start bp-end gen-map-filename\n";
		ShowEthnicities(connection);
	}

}
