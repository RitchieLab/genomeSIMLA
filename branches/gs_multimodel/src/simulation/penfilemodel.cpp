//
// C++ Implementation: penfilemodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "penfilemodel.h"
#include <sstream>
#include "utility/utility.h"

using namespace Utility;

namespace Simulation {

namespace StatusModel {

PenFileModel::PenFileModel() : modelSize(0), freqThreshold(0.0001) {
}


PenFileModel::~PenFileModel()
{
}


void PenFileModel::SetFilename(const char *filename) {
	this->filename=filename;
}


 
void PenFileModel::GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci) {
	os<<"<P><H3>Penetrance Based Model</H3>\n";
	if (modelSize != diseaseLoci.size()) 
		os<<"<B>Invalid Locus List Specified</B>";

	ReportPenetranceTable(os);
	
}

void PenFileModel::GenerateReport( ostream &os, uint headerWidth) {
	os<<setw(headerWidth)<<right<<""<<"#"<<modelID<<" (";
	os<<filename<<")";
	
	os<<endl;
}

void PenFileModel::VerifyAlleleFreq(char allele, float observedFreq) {
	float freq = alleleFreq[allele];
	stringstream errMsg;

	if (freq < observedFreq - freqThreshold || freq > observedFreq + freqThreshold) {

		errMsg<<allele<<"\t"<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<freq<<"\t";
		errMsg<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<observedFreq<<"\t";
		errMsg<<setiosflags(ios::fixed | ios::showpoint)<<setprecision(4)<<freqThreshold<<"\n";

		errMsg<<"An inconsistency has been identified with the allele frequencies in the pool and those specified\n";
		errMsg<<"in the penetrance configuration file. There are two reasons this might occur:\n";
		errMsg<<"\t1) The major/minor alleles are flipped. Capitol letters represent the first allele, regardless of\n";
		errMsg<<"\t   it's frequency. \n";
		errMsg<<"\t2) The pool was too small to support a variance of "<<freqThreshold<<" if you are drawing datasets\n";
		errMsg<<"\t   from generation 0. If this is the case, try increasing the initial size of the pool.\n";
		errMsg<<"\t3) The penetrance used is completely inappropriate for the specified loci. Are you sure you specified\n";
		errMsg<<"\t   the correct loci in the genomesim configuration?\n\n";
		errMsg<<"Below are the allele frequencies associated with the model loci:\n";
		
		char geno = 'A';
		char littleGeno = 'a';
		errMsg<<"\t"<<setw(10)<<"Locus"<<setw(10)<<"Obs. Freq."<<setw(10)<<"Specified Frequencies"<<"\n";
		for (uint i=0; i<loci.size(); i++) {
			errMsg<<"\t"<<"  Freq("<<loci[i].label<<","<<(char)(geno+i)<<") "<<setw(10)<<loci[i].alFreq1<<setw(10)<<(alleleFreq[(char)(geno+i)])<<"\n";
			errMsg<<"\t"<<"  Freq("<<loci[i].label<<","<<(char)(littleGeno+i)<<") "<<setw(10)<<loci[i].alFreq2<<setw(10)<<alleleFreq[(char)(littleGeno+i)]<<"\n";
		}
		throw Utility::Exception::General(errMsg.str().c_str());
	}
}

void PenFileModel::Refresh(PoolManager *pools) {
	size_t lociCount = loci.size();
	
	for (size_t i=0; i<lociCount; i++) {
		DiseaseLocus &locus = loci[i];
		pools->GetAlleleFrequency(locus.chromosome, locus.locusIdx, locus.alFreq1, locus.alFreq2);
	}
}



bool PenFileModel::Init(istream& s, PoolManager *pools) { 
	string labelType;
	bool indexBased = false;

	s>>labelType;

	filename = ParseFilename(s);

	if (labelType == "INDEX") {
		indexBased = true;
		DiseaseLocus l;
		bool doContinue = true;
		while (doContinue) {
			l = GetLocus(s);
			doContinue = l.ValidIdx();
			if (doContinue) {
				float a1, a2;
				if (pools->GetAlleleFrequency( l.chromosome, l.locusIdx, a1, a2))
					AddDiseaseLoci(l.label.c_str(), l.chromosome, l.locusIdx, a1);
				else {
					cerr<<"Bad locus specification  Chromosome: "<<l.chromosome<<" Locus: "<<l.locusIdx<<"\n";
					cerr<<"Be sure you are specifying chromosome and index\n";	
					return false;
				}
			}
		}
	}
	else if (labelType != "LABEL") {
		assert(0);
		return false;
	}
	else {
		DiseaseLocus l;
		bool doContinue = true;
			
		string label;
		while (doContinue) {
			label="";
			s>>label;

			doContinue = label.length() > 0;
			if (doContinue) {
				Locus loc;
				if (pools->GetLocusReference(label.c_str(), loc)) 
					AddDiseaseLoci(loc.GetLabel().c_str(), loc.GetChromID(), loc.GetID(), loc.GetMinAlleleFreq());
				else
					cout<<"Unable to find SNP, "<<label<<", anywhere!\n";
			}
		}
	}
	modelSize = loci.size();
	SetModelSize( modelSize );
	
	return true; 
}



bool PenFileModel::ParseLine(const char *line, uint val) {
	bool validEntry = false;
	if (line[0] > 64 && line[0] < 123) {
		stringstream ss(line);
		string key;
		ss>>key;

		if (strcmp(key.c_str(), LABEL_DISEASELOCI) == 0)   {
			//validEntry = true;
			//SetupDiseaseLoci(line);
		}
		else if (strcmp(key.c_str(), LABEL_PENTABLE) == 0) {
			validEntry = true;
		}
		else if (strcmp(key.c_str(), LABEL_FREQ_THRESH) == 0) {
			freqThreshold = 0.0;
			ss>>freqThreshold;
			if (freqThreshold == 0.0) {
				cerr<<"Invalid value for FREQ_THRESH: "<<line<<"\n";
				cerr<<"Unable to load the model.\n";
				abort();
			}
			validEntry = true;
		}
		else if (strcmp(key.c_str(), LABEL_FREQ) == 0) {
			char allele;
			double freq = 0.0;
			ss>>allele>>freq;
			alleleFreq[allele] = freq;
			validEntry = true;
		}
		else {
			if (modelSize > 0) {
				float penetrance;
				ss>>penetrance;
				AddPenetrance(GetGenotypeIdx(key.c_str(), modelSize), penetrance);
				validEntry = true;
			}
			else
				cout<<"Attempting to add a Penetrance to a model of size, 0.\n";
		}
	}
	return validEntry;
}

void PenFileModel::Load() {
	Load(filename.c_str());
}
string PenFileModel::GetModelConfiguration() {
	stringstream ss;
	ss<<"DEFINE_MODEL PENTABLE INDEX \""<<filename<<"\" ";
	ModelLociArray::iterator itr = loci.begin();
	ModelLociArray::iterator end = loci.end();

	while (itr != end) {
		cout<<".";
		ss<<itr->chromosome+1<<" "<<itr->locusIdx+1<<"\t";
//		ss<<itr->label<<" ";
		itr++;
	}
	cout<<"Model: "<<ss.str()<<"\n";

	cout<<"Model size: "<<loci.size()<<"\n";
	return ss.str();
}

void PenFileModel::Load(const char *filename) {
	LineParser lp;
	this->filename = filename;

	alleleFreq.clear();
	
	uint validLines = lp.Parse(filename, this);
	
	if (penCount < 3 || validLines < 1 + penCount) {
		cerr<<"Problems were encountered attempting to load the model file: "<<filename<<".\nPlease check that it is complete and ordered properly.\n";
		abort();
	}
	
	if (alleleFreq.size()  != 2 * loci.size() ){
		cout<<"As of genomesimla 2.0.4, all penetrance tables require the frequencies of each allele in letter format:\n";
		cout<<"\tFor instance\n\t\tFREQ A 0.8\n";
		throw Utility::Exception::General("Icorrectly configured penetrance model file. Please make sure that both alleles indicate the appropriate frequency (i.e. FREQ A 0.2) on separate lines.");
	}
		
	
	cout<<"\nDisease model in use: "<<filename<<"\n";
	int headerWidth = 45;
	
	cout<<setw(headerWidth-15)<<""<<loci[0].chromosome+1<<":"<<loci[0].locusIdx+1;
	for (uint i=1; i<modelSize; i++) 
		cout<<"x"<<loci[i].chromosome+1<<":"<<loci[i].locusIdx+1<<"\n";
	cout<<endl;
		
	cout<<"          --Penetrance table: \n";
	char geno = 'A';
	char littleGeno = 'a';
	for (uint i=0; i<loci.size(); i++) {
		VerifyAlleleFreq(geno+i, loci[i].alFreq1);
		VerifyAlleleFreq(littleGeno+i, loci[i].alFreq2);

		cout<<"\tFreq("<<(char)(geno+i)<<"): "<<loci[i].alFreq1<<"\n";
	}
	vector<string> pentable;
	BuildGenotypeLabels(pentable, modelSize);

	//We want to avoid meaningless tables like all zeros or 1.0s.
	uint fixedCells = 0;

	for (uint i=0; i<penCount; i++)	{
		cout<<"              "<<pentable[i]<<" : "<<penList[i]<<"\n";
		if	(penList[i] < 0.00000001 || penList[i] > 0.99999999)
			fixedCells++;
	}
	
	if (fixedCells == penCount) {
		throw Utility::Exception::General("The Penetrance table specified seems broken (it seems impossible to generate both Cases and Controls). Please verify that the table is correct. If it is, please contact us at genomesimla@chgr.mc.vanderbilt.edu to report a problem.");
	}

		

	cout<<endl;

}


}

}
