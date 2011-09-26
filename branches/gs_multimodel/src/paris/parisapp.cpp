/* 
 * File:   parisapp.cpp
 * Author: torstees
 * 
 * Created on January 5, 2010, 1:29 PM
 */

#include "parisapp.h"
#include "utility/strings.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "strings.h"
#include <algorithm>
#include "parislogs.h"
#include "bin.h"
#include "magicdb.h"

using namespace std;
using namespace soci;
using namespace Utility;

namespace Paris {

ParisResults ParisResults::db;
bool ParisResults::resultsDB = false;


ParisApp::ParisApp() {
}

ParisApp::ParisApp(const ParisApp& orig) {
}

ParisApp::~ParisApp() {
	std::map<std::string, Chromosome*>::iterator itr = chromosomes.begin();
	std::map<std::string, Chromosome*>::iterator end = chromosomes.end();

	while (itr != end) {
		delete itr++->second;
	}


}


void ParisApp::InitKnowledge(const char *dbFilename) {
	//Test the DB connection
	if (!Utility::FileExists(dbFilename)) {
		cerr<<"The database, "<<dbFilename<<", could not be found. Unable to continue.\n";
		exit(1);
	}
	try {
		this->dbFilename				= dbFilename;
		string cnxParam = "dbname="+string(dbFilename)+" timeout=10";
		sociDB.open(soci::sqlite3, cnxParam.c_str());


		string dbSnp, ensembl, hapmap, build, variations;
		sociDB<<"SELECT version FROM versions WHERE element='ncbi'", into(dbSnp);
		sociDB<<"SELECT version FROM versions WHERE element='ensembl'", into(ensembl);
		sociDB<<"SELECT version FROM versions WHERE element='hapmap'", into(hapmap);
		sociDB<<"SELECT version FROM versions WHERE element='variations'", into(variations);
		sociDB<<"SELECT version FROM versions WHERE element='build'", into(build);
		//this->varVersion				= atoi(variations.c_str());

		cout<<"\n------------------------- Dependency Versions ----------\n";
		cout<<setw(35)<<right<<"dbSNP: "<<dbSnp<<"\n";
		cout<<setw(35)<<right<<"Ensembl: "<<ensembl<<"\n";
		cout<<setw(35)<<right<<"Hap Map LD: "<<hapmap<<"\n";
		cout<<setw(35)<<right<<"Variations: "<<variations<<"\n";
		cout<<setw(35)<<right<<"Build: "<<build<<"\n";

	} catch (soci::soci_error const &e) {
		cerr<<"Problems were encountered trying to open the database, "<<dbFilename<<". Error: "<<e.what()<<"\n";
	}

}



uint ParisApp::InitBins(uint idealBinSize) {
	std::set<Feature*, SortByFeatureSize> binpool;
	std::set<Feature*> singleFeaturePool;

	//if (FileLogs::WriteBinDetails)
	FileLogs::logger.binDetails<<"Chrom,Gene ID,Start,Stop,RS ID,pvalue\n";

	std::map<std::string, Chromosome*>::iterator itr = chromosomes.begin();
	std::map<std::string, Chromosome*>::iterator end = chromosomes.end();

	while (itr != end) {
		itr->second->InitBins(binpool, singleFeaturePool, FileLogs::logger.binDetails);
		itr++;
	}
	uint activeBinCount = binpool.size();

	std::set<Feature*>::iterator sfItr = singleFeaturePool.begin();
	std::set<Feature*>::iterator sfEnd = singleFeaturePool.end();
	vector<Feature*> curBin;
	curBin.reserve(activeBinCount/4);
	uint spf = 0;

	while (sfItr != sfEnd) {
		Feature *f = *sfItr++;
		f->BinIndex(0);
		curBin.push_back(f);
		if (f->_begin == f->_end)
			spf++;

	}

	//We want to mix the bins up, so that our bins aren't grouped the
	//same for each random seed
	vector<Feature*> binShuffler;
	binShuffler.insert(binShuffler.begin(), binpool.begin(), binpool.end());
	std::random_shuffle(binShuffler.begin(), binShuffler.end(), Bin::rnd);
	binpool.clear();
	binpool.insert(binShuffler.begin(), binShuffler.end());

	cerr<<"Single Point Features: "<<spf<<"\n";

	std::set<Feature*, SortByFeatureSize>::iterator fItr = binpool.begin();
	std::set<Feature*, SortByFeatureSize>::iterator fEnd = binpool.end();

	bins.clear();
	bins[0] = curBin;

	//By doing this, we make sure it's the smallest deviation from the desired size
	uint binCount = (uint)(((float)(activeBinCount))/(float)idealBinSize+0.5);
	if (binCount == 0)
		binCount = 1;
	uint binSize = (uint)(((float)(activeBinCount)+0.5) / binCount);
	int sparePoints = (activeBinCount)%binCount;



	cerr<<setw(10)<<"Bin"<<setw(10)<<"Feature Count"<<setw(12)<<"Avg Bin Size"<<"\n";
	cerr<<setw(10)<<"0"<<setw(10)<<curBin.size()<<setw(12)<<1<<"\n";

	for (uint i=1; i<=binCount; i++) {
		uint minFeatureSize = (uint)-1;
		uint maxFeatureSize = 0;
		curBin.clear();
		uint localBinSize = binSize;
		if (sparePoints-- > 0)
			localBinSize++;
		curBin.reserve(localBinSize);

		double sumFeatureSize = 0.0;

		//stringstream blockBuffer;
		//int superBlock = 0;

		for (uint n=0; n<localBinSize && fItr != fEnd; n++) {
			Feature *f = *fItr++;
			if (f->FeatureSize() < minFeatureSize)
				minFeatureSize = f->FeatureSize();
			maxFeatureSize = f->FeatureSize();
			sumFeatureSize += (double)maxFeatureSize;


			/*
			 if (f->FeatureSize() > 200) {
				Chromosome *chr = chromosomes[f->_chromosome];
				newFile<<"<TABLE  CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES><TR><TH>Chromosome</TH><TH>Feature Beginning<TH><TH>Feature End</TH><TH>Feature Size</TH><TH>Link To SNPs</TH></TR>\n";

				newFile<<"\t<TR><TD>"<<f->_chromosome<<"</TD><TD>"<<f->_begin<<"</TD><TD>"<<f->_end<<"</TD><TD>"<<f->FeatureSize()<<"</TD><TD><A HREF=supersized-blocks.html#block_"<<superBlock<<"> Block "<<superBlock<<"</A></TD></TR>\n";
				blockBuffer<<"\n<A name=\"block_"<<superBlock<<"\"><H2>Block "<<superBlock<<"</H2></A><TABLE CELLSPACING=1 CELLPADDING=3 BORDER=1 RULES=ALL FRAME=HSIDES><TR><TH>Chromosome</TH><TH>RS ID</TH><TH>Position</TH><TH>P Value</TH><TH>Link</TH></TR>\n";
				superBlock++;
				map<uint, float> pvalues = f->GetPValues();
				map<uint, float>::iterator itr = pvalues.begin();
				map<uint, float>::iterator end = pvalues.end();
				while (itr != end) {
					blockBuffer<<"<TR><TD>"<<f->_chromosome<<"</TD><TD>"<<chr->GetSNP(itr->first)<<"</TD><TD>"<<itr->first<<"</TD><TD>"<<itr->second<<"</TD><TD><A HREF='http://uswest.ensembl.org/Homo_sapiens/Variation/Summary?source=dbSNP;v=rs"<<chr->GetSNP(itr->first)<<"'>rs"<<chr->GetSNP(itr->first)<<"</A></TD></TR>\n";
					itr++;
				}
				blockBuffer<<"</TABLE>\n";
			}*/
			
			f->BinIndex(i);
			assert(f->BinIndex() == i);
			curBin.push_back(f);
		}

		//newFile<<blockBuffer.str()<<"</HTML>\n";
		bins[i] = curBin;
		cerr<<setw(10)<<i<<setw(10)<<curBin.size()<<setw(12)<<(sumFeatureSize/(double)curBin.size())<<" ["<<minFeatureSize<<" - "<<maxFeatureSize<<"]"<<"\n";
	}


	double sumFeatureSize = 0.0;
	while (fItr != fEnd) {
		Feature *f = *fItr++;
		sumFeatureSize += (double)f->FeatureSize();
		f->BinIndex(binCount);
		assert(f->BinIndex() == binCount);
		bins[binCount].push_back(f);
	}
	if (bins[binCount].size() > binSize)
		cerr<<setw(10)<<binCount<<setw(10)<<bins[binCount].size()<<setw(12)<<(sumFeatureSize/(double)bins[binCount].size())<<"\n";


	itr = chromosomes.begin();

	while (itr != end) {
		itr->second->VerifyBins();
		itr++;
	}

	return bins.size();
}

void ParisApp::InitData(vector<ParisApp::DataPoint>& data, uint binSize) {
	vector<ParisApp::DataPoint>::iterator itr = data.begin();
	vector<ParisApp::DataPoint>::iterator end = data.end();
	cerr<<"Populating Data\n";

	uint observedPoints = 0,
		 ignoredPoints = 0;
	while (itr != end) {
		chromosomes[itr->chrom]->AddValue(itr->rsID, itr->pvalue);
		if (itr->pvalue > 0.0 || !Feature::IgnorePValueOfZero) {
			observedPoints++;
			//chromosomes[itr->chrom]->AddValue(itr->rsID, itr->pvalue);
		}
		else
			ignoredPoints++;
		itr++;
	}
	cerr<<"Observed :"<<observedPoints<<"\tIgnored Points: "<<ignoredPoints<<"\n";
	cerr<<"Binning\n";
	//Perform Binning
	InitBins(binSize);
}

std::set<Pathway*> ParisApp::GetPathways(std::multimap<uint, uint>& groups) {
	std::multimap<uint, uint>::iterator itr = groups.begin();
	std::multimap<uint, uint>::iterator end = groups.end();
	std::set<Pathway*> pathways;
	while (itr != end) {
		if (knowledge.find(itr->first) != knowledge.end())
			pathways.insert(knowledge[itr->first]->GetPathway(itr->second));
		itr++;
	}
	return pathways;
}

void ParisApp::InvestigatePathway(uint permutationCount, float dataSig, float pathSig, const char *pathwayName, bool showAllPathways) {

	if (string(pathwayName).length() > 0) {
		std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
		std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();
		cerr<<"Investigating Pathway "<<pathwayName<<":\n";
		Pathway *pathway = NULL;
		KnowledgeBase *kb = NULL;

		map<uint, Analyzer::Result> resultBuffer;

		while (kbItr != kbEnd && pathway == NULL) {
			kb = kbItr->second;
			pathway = kb->GetPathway(pathwayName);
			//kbItr->second->InvestigatePathways(bins, permutationCount, significance, scores, pathwayName);
			//kbItr->second->ReportResults(scores, report, true);
			kbItr++;
		}

		if (pathway) {
			ofstream snpReport(FileLogs::logger.GetFilename((pathway->Name() + "-snp-report").c_str(), "csv").c_str());
			pathway->GenerateSnpReport(snpReport,chromosomes, true);
			snpReport.close();
			string filename = FileLogs::logger.GetFilename(pathway->Name().c_str(), "html");
			ofstream file(filename.c_str());
			file<<
	"<HTML><HEAD><TITLE>PARIS Pathway Breakdown "<<reportName<<" "<<pathwayName<<"</TITLE></HEAD>\n"<<
	"<BODY LINK=#4F2412 VLINK=#4F2412>\n";
	


			std::set<Gene*> genes = pathway->GetGenes();
			std::set<Gene*>::iterator gItr = genes.begin();
			std::set<Gene*>::iterator gEnd = genes.end();
			Analyzer analyzer;
			Analyzer::Result pathwayResult;
			if (resultBuffer.find(pathway->GroupID()) == resultBuffer.end()) {
				pathwayResult = analyzer.Analyze(pathway->GroupID(), genes, bins, permutationCount, false);
				resultBuffer[pathway->GroupID()] = pathwayResult;
			}
			else
				pathwayResult = resultBuffer[pathway->GroupID()];
			file<<"<H2>Pathway Investigation : ";
			if (reportName.length() > 0)
				file<<reportName<<" - ";
			if (kb->GetID() == 2)
				file<<"<A HREF='http://www.genome.jp/dbget-bin/www_bget?"<<pathwayName<<"' TEXT=A3C586>"<<pathwayName<<"</A> </H2><P>"<<pathway->Description()<<"<P>\n";
			else
				cerr<<"ID "<<kb->GetID()<<" </H2><P>"<<kb->GetName()<<"\n";
			file<<"<P><TABLE border=1 cellspacing=0>"
				 <<"<TR bgcolor=#CCCCCC border=1 styles='color:#A3C586;'><TH>Feature Count</TH>"
				 <<"<TH>Simple Features</TH><TH>Simple Features*</TH>"
				 <<"<TH>Complex Features</TH><TH>Complex Features*</TH>"
				 <<"<TH>P Value n/"<<permutationCount<<"</TH></TR>\n";
			file<<"\t<TD>"<<pathwayResult.complexFeatures+pathwayResult.simpleFeatures+pathwayResult.sigSimple+pathwayResult.sigComplex
				 <<"</TD><TD>"<<pathwayResult.simpleFeatures<<"</TD><TD>"<<pathwayResult.sigSimple
				 <<"</TD><TD>"<<pathwayResult.complexFeatures<<"</TD><TD>"<<pathwayResult.sigComplex<<"</TD><TD>"<<pathwayResult.GetPValue()<<"</TD></TR></TABLE>\n";
			file<<"\t<H2>Gene Breakdown for Pathway "<<pathwayName<<"</H2>\n<CENTER>\n";
			file<<"\t<P><TABLE border=1 cellspacing=0><TR bgcolor=#CCCCCC border=1><TH>Gene Name</TH><TH>Ensembl ID</TH><TH>Feature Count</TH><TH>Simple Feature</TH>"<<
				 "<TH>Simple Feature*</TH><TH>Complex Features</TH><TH>Complex Features*</TH><TH>Gene PValue n/"<<permutationCount<<"</TH><TH>Associated Pathways</TH></TR>\n";

			while (gItr != gEnd) {
				Gene *gene = *gItr;
				std::set<Gene*> singleGene;
				singleGene.insert(gene);

				Analyzer::Result result = analyzer.Analyze(gene->id, singleGene, bins, permutationCount, false);
				if (result.PValue() <= pathSig)
					file<<"\t<TR bgcolor=#E9E0DB border=1>";
				else
					file<<"\t<TR border=1>";
				file<<"<TD>"<<gene->Alias()<<"</TD><TD><A HREF='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g="<<gene->EnsemblID()<<"'>"<<gene->EnsemblID()<<"</A></TD><TD>"<<gene->FeatureCount()
					<<"</TD><TD>"<<result.simpleFeatures<<"</TD><TD>"<<result.sigSimple<<"</TD>"
					<<"<TD>"<<result.complexFeatures<<"</TD><TD>"<<result.sigComplex<<"</TD>"
					<<"<TD>"<<result.GetPValue()<<"</TD><TD>";
					std::multimap<uint, uint> groups = gene->GetGroups();
					std::set<Pathway*> pathways = GetPathways(groups);
					std::set<Pathway*>::iterator itr = pathways.begin();
					std::set<Pathway*>::iterator end = pathways.end();
					int count = 0;


					while (itr != end) {
						if (count++ > 0)
							file<<" ";
						std::set<Gene*> pwGenes = (*itr)->GetGenes();

						Analyzer::Result pwResult;
						if (resultBuffer.find((*itr)->GroupID()) == resultBuffer.end()) {
							pwResult = analyzer.Analyze((*itr)->GroupID(), pwGenes, bins, permutationCount, false);
							resultBuffer[(*itr)->GroupID()] = pwResult;
						}
						else
							pwResult = resultBuffer[(*itr)->GroupID()];
						if (pathway->GroupID() != (*itr)->GroupID()) {
							if (pwResult.PValue() < pathSig)
								file<<"<B>"<<(*itr)->Name()<<"</B>";
							else if (showAllPathways)
								file<<(*itr)->Name();
						}
						itr++;
					}
					file<<"</TD></TR>\n";

				gItr++;
			}
			file<<"</TABLE><P><I>*pvalues of "<<dataSig<<" or less.<P>Shaded rows indicate genes whose pvalues are "<<pathSig<<" or less.</BODY></HTML>\n";
			cerr<<pathway->Name()<<" : "<<filename<<"\n";
		}
		else
			cerr<<"Unable to find pathway: "<<pathwayName<<". Is the pathway name spelled correctly, and is the knowledge base active?\n";
	}
	//Evaluate each group for statistical value (including the ptests)
}

void ParisApp::WriteKbVitalsReport(const char *vitalsReport) {
	ofstream file(vitalsReport);
	file<<"kb\tpathway-name\tpathway-id\tgene-size\tgene-size-sd\tgene-contents\tgene-contents-sd\tavg-bl-size\tavg-bl-size-sd\n";
	std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
	std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();
	while (kbItr != kbEnd) {
		kbItr->second->WritePathwayStatistics(file);
		kbItr++;
	}
}

void ParisApp::WriteBinReport(const char* name) {
	std::map<std::string, Chromosome*>::iterator itr = chromosomes.begin();
	std::map<std::string, Chromosome*>::iterator end = chromosomes.end();

	ofstream file(name);
	while (itr != end) {
		itr->second->WriteBinReport(file);
		itr++;
	}
}

map<std::string, uint> ParisApp::AnalyzeNegativeControl(uint permutationCount, float datasetSignificance, float pathSignificance) {
	std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
	std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();

	map<std::string, uint> sigTests;
	float sigSum = 0.0;

	if (kbItr != kbEnd) {
		if (reportName.length() > 0)
			FileLogs::logger.finalReport<<reportName<<" ";
		kbItr->second->ResultsHeader(FileLogs::logger.finalReport);
		while (kbItr != kbEnd) {
			std::multiset<Analyzer::Result> scores;
			if (FileLogs::WriteDetailedFeatureReport)
				kbItr->second->DetailedReport(chromosomes, FileLogs::logger.detailedFeatureReport);
			uint sigs = kbItr->second->AnalyzeNegativeControl(bins, permutationCount, pathSignificance, scores);
			sigSum+=(float)sigs;
			sigTests[kbItr->second->GetName()] = sigs;
			kbItr->second->ReportResults(scores, FileLogs::logger.finalReport);
			cerr<<"--"<<kbItr->second->GetName()<<"\tSignificant Tests:"<<sigs<<"\n";
			kbItr++;

		}
		knowledge.begin()->second->ResultsFooter(FileLogs::logger.finalReport, datasetSignificance, pathSignificance);
		cerr<<"\nAverage Significant Count: "<<(sigSum/(float)(sigTests.size()))<<"\n";


	}
	return sigTests;
	//Evaluate each group for statistical value (including the ptests)
}

void ParisApp::RunAnalysis(uint permutationCount, float datasetSignificance, float pathSignificance) {
	std::map<uint, KnowledgeBase*>::iterator kbItr = knowledge.begin();
	std::map<uint, KnowledgeBase*>::iterator kbEnd = knowledge.end();
	if (kbItr != kbEnd) {
		if (reportName.length() > 0)
			FileLogs::logger.finalReport<<reportName<<" ";
		kbItr->second->ResultsHeader(FileLogs::logger.finalReport);
		while (kbItr != kbEnd) {
			std::multiset<Analyzer::Result> scores;
			if (FileLogs::WriteDetailedFeatureReport)
				kbItr->second->DetailedReport(chromosomes, FileLogs::logger.detailedFeatureReport);
			kbItr->second->RunAnalysis(bins, permutationCount, pathSignificance, scores);
			kbItr->second->ReportResults(scores, FileLogs::logger.finalReport);
			kbItr++;
		}
		knowledge.begin()->second->ResultsFooter(FileLogs::logger.finalReport, datasetSignificance, pathSignificance);
	}

	//Evaluate each group for statistical value (including the ptests)
}

void ParisApp::InitKB(const char* popID, uint geneExpansion, Utility::StringArray& groups) {
		if (ParisResults::resultsDB) {
			string resultsDB = reportPrefix + "-results.sqlite";
			ParisResults::db.Open(resultsDB.c_str());
			ParisResults::db.InitTable("knowledgebase",	"CREATE TABLE knowledgebase (kb_type INTEGER, kb_name VARCHAR(40))");
			ParisResults::db.InitTable("pathways",			"CREATE TABLE pathways (pathway_id INTEGER, kb_type INTEGER, pathway_name VARCHAR, pathway_description VARCHAR)");
			ParisResults::db.InitTable("genes",				"CREATE TABLE genes (gene_id INTEGER, ensembl_id VARCHAR, chrom VARCHAR, start INTEGER, end INTEGER)");
			ParisResults::db.InitTable("features",			"CREATE TABLE features (feature_id INTEGER, chrom STRING, start INTEGER, end INTEGER)");
			ParisResults::db.InitTable("feature_snps",	"CREATE TABLE feature_snps (feature_id INTEGER, rsid INTEGER, pos INTEGER, pvalue DECIMAL)");
			ParisResults::db.InitTable("gene_to_feature", "CREATE TABLE gene_to_feature (gene_id INTEGER, feature_id INTEGER)");
			ParisResults::db.InitTable("pathway_to_gene", "CREATE TABLE pathway_to_gene (pathway_id INTEGER, gene_id INTEGER)");
			ParisResults::db.InitTable("snps",				"CREATE TABLE snps (chrom STRING, rsid INTEGER, pos INTEGER)");
		}

		cerr<<"Initializing knowledge for population: "<<popID<<"\n";
		std::map<string, Chromosome*>::iterator itr = chromosomes.begin();
		std::map<string, Chromosome*>::iterator end = chromosomes.end();
		uint featureID = 0;

		//Load the alias lookup table
		rowset<row> rs = (sociDB.prepare <<"SELECT DISTINCT alias, gene_id FROM region_alias NATURAL JOIN regions WHERE region_alias_type_id=1300 GROUP BY gene_id");
		map<uint, string> aliasLookup;
		for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
			row const& row = *itr;
			string alias = row.get<string>(0);
			uint geneID = row.get<int>(1);
			aliasLookup[geneID] = alias;
		}

		while (itr != end) {
			//First step, is to load the features
			itr->second->LoadFeatures(sociDB, popID, featureID);
			//Next, we want to Load the genes, and then merge the features into the genes
			itr->second->LoadGenes(sociDB, 0, geneExpansion, aliasLookup);
			//itr->second->MergeFeaturesIntoGenes();
			itr++;
		}

		Utility::StringArray::iterator gitr = groups.begin();
		Utility::StringArray::iterator gend = groups.end();

		while (gitr != gend) {
			uint groupType=0;
			string groupName="";
	//		try {
				//Finally, we want to set up the group information and associated genes with the various groupings
				sociDB<<"SELECT group_type_id, group_type FROM group_type WHERE group_type_id=?", use(atoi(gitr->c_str())), into(groupType), into(groupName);
				if (ParisResults::resultsDB) 
					ParisResults::db.sociDB<<"INSERT INTO knowledgebase VALUES (:type, :name)", use(groupType), use(groupName);
				KnowledgeBase *kb=new KnowledgeBase(groupType, groupName.c_str());

				cerr<<"Loading Knowledgebase "<<groupType<<" "<<groupName<<"\n";
				kb->LoadKnowledge(sociDB, chromosomes, FileLogs::logger.kbReport);
				knowledge[groupType] = kb;
/*			} catch (soci_error const &e) {
				cerr<<"Unable to Load group data from database. Error: "<<e.what()<<"\n";
				exit(1);
			}
*/			gitr++;
		}
	
}

void ParisApp::ListGroupIDs() {
	rowset<row> rs = (sociDB.prepare <<"SELECT DISTINCT group_type_id, group_type, download_date FROM group_type");
	cerr<<"\n"<<setw(10)<<right<<"KB ID"<<setw(30)<<right<<"Name"<<setw(30)<<right<<"Date "<<"\n";
	for (rowset<row>::const_iterator itr = rs.begin(); itr != rs.end(); ++itr) {
		row const& row = *itr;
		uint groupID = row.get<int>(0);
		string name = row.get<string>(1);
		std::tm downloaddate = row.get<std::tm>(2);
		std::time_t time = std::mktime(&downloaddate);

		cerr<<setw(10)<<right<<groupID<<setw(30)<<right<<name<<setw(30)<<right<<std::ctime(&time);

	}

}

uint ParisApp::InitSNPs(std::multimap<string, uint>& allSNPs) {
	string filename;
	string path, basename, extension;
	Utility::SplitIntoComponents(dbFilename.c_str(), path, basename, extension);

	sociDB<<"SELECT version FROM versions WHERE element = 'variations'", into(filename);

	if (filename.length() == 0) {
		cerr<<"There is a problem with the database. It doesn't have a link to a variations file.";
		abort();
	}

	if (!Utility::FileExists(filename.c_str())) {
		filename = path + "/" + filename;
	}

	
	
	
	
	
	ifstream file(filename.c_str(), ios::binary);
	if (!file.good()) {
		cerr<<"Unable to open file, "<<filename<<". Aborting\n";
		exit(1);
	}
	uint count = 0;
	int varVersion = 0;
	file.read((char*)&varVersion, 4);
	
	while (file.good()) {
		char label[3];
		file.read(label, 2);

		if (!file.eof()) {
			label[2]='\0';

			string chrID = label;
			chrID.erase(remove_if(chrID.begin(), chrID.end(), ::isspace), chrID.end());


			int snpCount=0, maxPosition=0;
			file.read((char*)&snpCount, 4);
			file.read((char*)&maxPosition, 4);
			//cout<<"+++ "<<chrID<<" "<<snpCount<<" "<<maxPosition<<"\n";
			std::multimap<string, uint>::iterator chromItr = allSNPs.lower_bound(chrID);
			std::multimap<string, uint>::iterator chromEnd = allSNPs.upper_bound(chrID);

			std::set<uint> snps;
			while (chromItr != chromEnd)  {
				snps.insert(chromItr->second);
				chromItr++;
			}

			Chromosome *newChrom = new Chromosome(label);

	//		offset+= maxPosition;
	//		posLookup[offset] = newChrom;
			if (file.good()) {
				cerr<<".";cerr.flush();
				for (int i=0; i<snpCount; i++) {
					uint rs=0, pos=0, role=0;
					file.read((char*)&rs, 4);
					file.read((char*)&pos, 4);
					file.read((char*)&role, 1);
					if (snps.size() == 0 || (snps.find(rs) != snps.end())) {
						//EST-RS>0 if (rs > 0) {
							newChrom->AddSNP(rs, pos);
							//cout <<" ++++ "<<rs<<" "<<pos<<" "<<role<<"\n";
							count++;
						//}
					}
				}
			}
			chromosomes[chrID] = newChrom;
		}
	}

	return count;
}

void ParisApp::SetReportPrefix(const char *prefix) {
	reportPrefix = prefix;
}

void ParisApp::ReportName(const char *name) {
	reportName = name;
}
std::string ParisApp::ReportName() {
	return reportName;
}

}
