/* 
 * File:   ldspline.cpp
 * Author: torstees
 * 
 * Created on August 30, 2010, 12:16 PM
 */

#include "ldspline.h"
#include <iostream>
#include <sstream>
#include <iomanip>


LdSplineApp::LdSplineApp()  {}
LdSplineApp::LdSplineApp(const LdSplineApp::LdSplineApp& orig) {}
LdSplineApp::~LdSplineApp() {}

void LdSplineApp::LoadHeadersFromHapmap(const char *chrom, const char* filename) {
	LocusLookup loci(&file, chrom);
	filenames[chrom] = filename;
	loci.LoadHeadersFromHapmap(filename);
	this->loci[chrom] = loci;
}
void LdSplineApp::ReportBounds(const char *chrom, int position, float value, const char *type, std::ostream& os) {
	if (loci.find(chrom) != loci.end()) {
		std::pair<int, int> bounds(position, position);

		if (strcmp(type, "dp") == 0 || strcmp(type, "DP") == 0)
			bounds = loci[chrom].GetLdSplineBoundsDP(position, value);
		else if (strcmp(type, "rs") == 0 || strcmp(type, "RS") == 0)
			bounds = loci[chrom].GetLdSplineBoundsRS(position, value);
		else {
			std::cerr<<"Unknown LD type: "<<type<<"\n";
			abort();
		}

		os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\trs"
			 <<loci[chrom].GetPosToRS(bounds.first)<<"\t"<<bounds.first<<"\trs"
			 <<loci[chrom].GetPosToRS(bounds.second)<<"\t"<<bounds.second<<"\t"<<value<<"\n";
	}
	else {
		os<<chrom<<"\t ?? \t"<<position<<"\t ??\t"<<position<<"\t ??\t"<<position<<"\n";
	}
}
void LdSplineApp::RunReport(const char *chrom, int position, float value, const char *type, std::ostream& os) {
	if (loci.find(chrom) != loci.end()) {
		std::map<int, float> spline;

		if (strcmp(type, "dp") == 0 || strcmp(type, "DP") == 0)
			spline = loci[chrom].GetLdSplineDP(position, value);
		else if (strcmp(type, "rs") == 0 || strcmp(type, "RS") == 0)
			spline = loci[chrom].GetLdSplineRS(position, value);
		else {
			std::cerr<<"Unknown LD type: "<<type<<"\n";
			abort();
		}
		std::map<int, float>::iterator sitr = spline.begin();
		std::map<int, float>::iterator send = spline.end();

		while (sitr != send) {
			std::string direction  = "down";
			if (position < sitr->first)
				direction = "up";
			os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\trs"<<loci[chrom].GetPosToRS(sitr->first)<<"\t"<<sitr->first<<"\t"<< std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setprecision(4)<<sitr->second<<"\t"<<direction<<"\n";
			sitr++;
		}

		if (spline.size() == 0)
			if (loci[chrom].GetPosToRS(position) > 0)
				os<<chrom<<"\trs"<<loci[chrom].GetPosToRS(position)<<"\t"<<position<<"\t \t \n";
			else
				os<<chrom<<"\t  ??  \t"<<position<<"\t \t \n";
	}
	else {
		os<<chrom<<"\t ?? \t"<<position<<"\t"<<"??"<<"\t"<<" "<<"\n";
	}
}

std::vector<SnpSpline> LdSplineApp::GetLocusRange(const char *chrom, int start, int stop) {
	std::vector<SnpSpline> splines;
	if (loci.find(chrom) != loci.end()) {
		splines = loci[chrom].GetLocusRange(start, stop);
	}
	return splines;
}

void LdSplineApp::Summarize(std::ostream& os, const char *chrom) {
	std::map<std::string, LocusLookup>::iterator litr = loci.begin();
	std::map<std::string, LocusLookup>::iterator lend = loci.end();
	while (litr != lend) {
		if (litr->first == chrom || litr->first == "ALL")
			litr->second.Summarize(os);
		litr++;
	}
}

void LdSplineApp::RunReport(std::ostream& os) {
	std::map<std::string, LocusLookup>::iterator litr = loci.begin();
	std::map<std::string, LocusLookup>::iterator lend = loci.end();

	os<<"Report\n";
	while (litr != lend) {
		LocusLookup::PosToRS rsLookup = litr->second.GetPosToRS();
		LocusLookup::PosToRS::iterator itr = rsLookup.begin();
		LocusLookup::PosToRS::iterator end = rsLookup.end();

		while (itr != end) {
			std::map<int, float> spline =  litr->second.GetLdSplineDP(itr->first, 0.90);
			std::stringstream st100;
			st100<<" DP 0.90: ";
			std::map<int, float>::iterator sitr = spline.begin();
			std::map<int, float>::iterator send = spline.end();
			while (sitr != send) {
				st100<<"\t"<<sitr->first<<" ("<<sitr->second<<")";
				sitr++;
			}

			spline =  litr->second.GetLdSplineDP(itr->first, 0.01);
			std::stringstream st75;
			st75<<" DP 0.01: ";
			sitr = spline.begin();
			send = spline.end();
			while (sitr != send) {
				st75<<"\t"<<sitr->first<<" ("<<sitr->second<<")";
				sitr++;
			}

			os<<itr->first<<"\t"<<itr->second<<"\n\t"<<st75.str()<<"\n\t"<<st100.str()<<"\n";

			//std::cerr<<itr->first<<"\t"<<itr->second<<"\t"<<loci.GetLdSplineDP(itr->first, 1.0).size()<<"\t"<<loci.GetLdSplineDP(itr->first, 0.9).size()<<"\n";
			itr++;
		}
		litr++;
	}
}

void LdSplineApp::SaveToCopyBinary(const char *newFilename) {
	file.flush();
	std::fstream newFile(newFilename,  std::ios::out|std::ios::in|std::ios::binary|std::ios::trunc);

	if (newFile.fail()) {
		std::cerr<<"Unable to open file, "<<newFilename<<"\n";
		abort();
	}

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();

	int count = loci.size();
//std::cerr<<"Chromosome Count: "<<count<<"\n";
	newFile.write((char*)&count, 4);
	int64_t offset = 4;

	while (itr != end)
		offset += (8 + itr++->second.LocusCount() * 16);

	// At this point, we should have the right offset for the first chromosome's splines
	itr = loci.begin();
	while (itr != end) {
		//itr->second.LoadBinaryHeader();
		offset = itr->second.DumpBinaryHeader(&newFile, offset);
		itr++;
	}

	itr = loci.begin();
	while (itr != end) {
		std::string filename = filenames[itr->first];
		//itr->second.LoadBinary();
		itr->second.DumpBinary(&newFile);
		itr++;
	}

	newFile.flush();
}

void LdSplineApp::SaveToBinary(const char *filename) {
	std::cerr<<"Opening file: "<<filename<<"\n";
	file.close();
	file.open(filename, std::ios::out|std::ios::in|std::ios::binary|std::ios::trunc);

	if (file.fail()) {
		std::cerr<<"Unable to open file, "<<filename<<"\n";
		abort();
	}
	//std::ofstream bin(filename, std::ios::binary);

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();

	int count = loci.size();
//std::cerr<<"Chromosome Count: "<<count<<"\n";
	file.write((char*)&count, 4);

	int64_t offset = 4;

	while (itr != end)
		offset += (8 + itr++->second.LocusCount() * 16);

	// At this point, we should have the right offset for the first chromosome's splines
	itr = loci.begin();
	while (itr != end) {
		offset = itr->second.DumpBinaryHeader(offset);
		itr++;
	}

	itr = loci.begin();
	while (itr != end) {
		std::string filename = filenames[itr->first];
		itr->second.LoadLdFromHapmap(filename.c_str());
		itr->second.DumpBinary();
		itr->second.Release();
		itr++;
	}

	file.flush();

}

void LdSplineApp::LoadFromBinary(const char *filename) {
	//std::ifstream bin(filename, std::ios::binary);
	file.open(filename, std::ios::in|std::ios::out|std::ios::binary);
	int count = 0;
	file.read((char*)&count, 4);

//std::cerr<<"Chromosome Count: "<<count<<"\n";
	for (int i=0; i<count; i++) {
		LocusLookup chromosome(&file);
		chromosome.LoadBinaryHeader();
		loci[chromosome.Chromosome().c_str()] = chromosome;
	}

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();
	while (itr != end) 
		itr++->second.LoadBinary();
}


void LdSplineApp::OpenBinary(const char *filename, bool loadFullHeaders) {
	file.open(filename, std::ios::in|std::ios::out|std::ios::binary);
	int count = 0;
	file.read((char*)&count, 4);

//std::cerr<<"Chromosome Count: "<<count<<"\n";
	for (int i=0; i<count; i++) {
		LocusLookup chromosome(&file);
		if (loadFullHeaders) 
			chromosome.LoadBinaryHeader();
		else
			chromosome.SkimBinaryHeader();
		//chromosome.LoadBinaryHeader();
		loci[chromosome.Chromosome().c_str()] = chromosome;
	}
}

void LdSplineApp::ExportForLiftOver(const char *bimOrig) {
	std::ofstream file(bimOrig);

	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();
	while (itr != end) 
		itr++->second.WriteMapForLiftOver(file);
}

void LdSplineApp::ImportLiftOver(const char *loFilename, const char *loUnmapped) {
	std::map<std::string, std::set<int> > unmappedPositions;		///< Chr -> position [position...]

	//Create a bunch of sets...multimap makes it annoying to do searches for the resulting contents
	std::map<std::string, LocusLookup>::iterator itr = loci.begin();
	std::map<std::string, LocusLookup>::iterator end = loci.end();


	std::ifstream infile(loUnmapped);
	char line[4096];
	while (infile.good() && !infile.eof()) {
		infile.getline(line, 4096);
		if (line[0]!= '#') {
			std::stringstream ss(line);
			std::string chr;
			int pos = 0;

			ss>>chr>>pos;

			if (pos > 0) {
				std::cerr<<pos<<" - unmapped\n";
				chr.erase(0,3);
				unmappedPositions[chr].insert(pos);
				//unmappedPositions.insert(std::pair<std::string, int>(chr, pos));
			}
		}
	}

	infile.close();
	infile.open(loFilename);

	itr = loci.begin();
	while (itr != end) {
		itr++->second.LoadMapFromLiftOver(infile, unmappedPositions[itr->first]);
	}
	
}

std::string LdSplineApp::ConvertHaploLD(const char *filename) {
	std::ifstream haplo(filename);

	while (haplo.good() && !haplo.eof()) {
		std::string chrom = "", filename = "";
		haplo>>chrom>>filename;

		if (chrom != "")
			LoadHeadersFromHapmap(chrom.c_str(), filename.c_str());
	}


	std::string newFilename = std::string(Utility::ExtractBaseFilename(filename)) + ".ldspline";
	SaveToBinary(newFilename.c_str());

	return newFilename;
}


int main(int argc, char **argv) {
	if (argc < 3) {
		std::cerr<<"Usage: ldspline [load/list/report/summarize/report-bounds] filename [dp/rs] [min-ld_value] [chrom] [position|filename]..\n";
		std::cerr<<"\tload         - initiates loading LD information in haploview format.\n";
		std::cerr<<"\tlist         - lists the spline results for a given set of positions. \n";
		std::cerr<<"\treport-bounds- lists spline boundaries for a give set of positions.\n";
		std::cerr<<"\tsummarize    - lists basic details about loci contained within the data store. Summarize requires a chromosome number of ALL\n";
		std::cerr<<"\t\tThe first parameter after list should be the file containing the LD data being queried.\n";
		std::cerr<<"\tdp/rs        - (list only, required) Indicates which type of statistic is being queried for (required only for list)\n";
		std::cerr<<"\tmin-ld_value - (list only, required) Indicates the threshold for the spline.\n";
		std::cerr<<"\tchrom        - (list only, required) Indicates which chromosome the positions are from.\n";
		std::cerr<<"\n* When running lists, users supply positions instead of filenames\n";
		return 1;
	}

	std::string cmd = argv[1];
	std::string dataFilename = Utility::StripExtension(argv[2]);
	std::string rawFilename = argv[2];

	if (cmd == "load") {
		for (int i=2; i<argc; i++) {
			LdSplineApp ldspline;
			std::string filename = ldspline.ConvertHaploLD(rawFilename.c_str());
			//ldspline.RunReport(std::cerr);
		}
	}

	else if (cmd == "range") {
		if (argc > 4) {
			LdSplineApp ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string chrom = argv[3];
			int first = atoi(argv[4]);
			int last  = atoi(argv[5]);

			ldspline.GetLocusRange(chrom.c_str(), first, last);
		}
		else
			std::cerr<<"You need more arguments for that!\n";
	}


	else if (cmd == "list") {
		if (argc > 6) {
			LdSplineApp ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string type = argv[3];
			float ldvalue = atof(argv[4]);
			std::string chrom = argv[5];
			std::cout<<"chrom\tprimary rs\tprimary pos\tsecondary rs\tsecondary pos\t"<<type<<"/tdirection\n";
			for (int i=6; i<argc; i++) {
				ldspline.RunReport(chrom.c_str(), atoi(argv[i]), ldvalue, type.c_str(), std::cout);
			}
		}
		else {
			std::cerr<<"-- Insufficient parameters for list\n";
		}
	}
	else if (cmd == "summarize") {
		LdSplineApp ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
		std::string chrom = "ALL";
		if (argc > 3)
			chrom = argv[3];
		ldspline.Summarize(std::cerr, chrom.c_str());
	}
	else if (cmd == "report-bounds") {
		if (argc > 6) {
			LdSplineApp ldspline;
			ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
			std::string type = argv[3];
			float ldvalue = atof(argv[4]);
			std::string chrom = argv[5];
			std::cout<<"chrom\trs\tpos\tlower_pos\tlower_rs\tupper_bound\tupper_rs\t"<<type<<"\n";
			for (int i=6; i<argc; i++) {
				ldspline.ReportBounds(chrom.c_str(), atoi(argv[i]), ldvalue, type.c_str(), std::cout);
			}
		}
		else {
			std::cerr<<"-- Insufficient parameters for list\n";
		}

	}
	else if (cmd == "report") {
		LdSplineApp ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
		ldspline.RunReport(std::cerr);
	}
	else if (cmd == "export-lomap") {
		LdSplineApp ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str());
		std::string loFilename = Utility::ExtractBaseFilename(argv[2]) + ".bim";
		ldspline.ExportForLiftOver(loFilename.c_str());
		std::cerr<<"Lift over output: "<<loFilename<<"\n";
	}
	else if (cmd == "import-lomap") {
		LdSplineApp ldspline;
		ldspline.OpenBinary(std::string(dataFilename + ".ldspline").c_str(), true);
		std::string loFilename = Utility::ExtractBaseFilename(argv[2]) + ".new";
		std::string loUnmapped = Utility::ExtractBaseFilename(argv[2]) + ".unmapped";
		std::cerr<<"Loading position data from "<<loFilename<<"\t"<<loUnmapped<<"\n";
		ldspline.ImportLiftOver(loFilename.c_str(), loUnmapped.c_str());
		ldspline.SaveToCopyBinary(std::string(dataFilename + "-b37.ldspline").c_str());
	}
	else if (cmd == "asdf") {
		for (int i=2; i<argc; i++) {
			LdSplineApp ldspline2;
			LdSplineApp ldspline;
			std::string filename = ldspline.ConvertHaploLD(argv[i]);
			ldspline2.LoadFromBinary(filename.c_str());
			ldspline.RunReport(std::cerr);
			ldspline2.RunReport(std::cerr);
		}
	}
	else {
		std::cerr<<"Unknown command: "<<argv[1]<<" Options include: load|list|summarize\n";
	}




	return 0;
}