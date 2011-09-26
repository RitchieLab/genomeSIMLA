// PenetranceModel.cpp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <math.h>
#include <iomanip>
#include "penetrancemodel.h"
#include "penfilemodel.h"
#include "gamodel.h"
#include "simlamodel.h"
#include "individual.h"
#include "utility/random.h"
#include "utility/exception.h"
#include "meantable.h"
#include "pedigreereferencesample.h"
//#include "defaults.h"

namespace Simulation {

namespace StatusModel {

using namespace std;
using namespace Utility;


const char *PenetranceModel::PenTable 	= "PENTABLE";
const char *PenetranceModel::Simla  	= "SIMLA";
const char *PenetranceModel::Simpen   	= "SIMPEN";
const char *PenetranceModel::Continuous	= "CONTINUOUS";
const char *PenetranceModel::Template	= "TEMPLATE";



uint PenetranceModel::GetMultiplier(uint genotype, uint position) {
	uint modelSize = loci.size();
	assert(genotype < 3);
	uint pwr=modelSize - position - 1;
	return (uint)((float)genotype * powf(3.0, (float)pwr));
}

void PenetranceModel::SetModelSize(size_t modelSize) {
	if (penList)
		delete[] penList;	
	penCount = (uint)powf(3.0, (modelSize));
	penList = new double[penCount];

}

void PenetranceModel::WritePenetranceFile(const char *filename, vector<Locus *>& diseaseLoci, float thresh) {
	ofstream file(filename);
	
	if (!file.is_open()) {
		throw Utility::Exception::FileNotWritable(filename);
	}
	
	//Need to work out a comment section	
	file<<"# The Threshold will help protect genomeSIMLA from using a penetrance\n";
	file<<"# table with data whose allele frequencies vary too greatly from the population\n";
	file<<"# for which the table was designed.\n";
	file<<"FREQ_THRESHOLD "<<thresh<<"\n";
	file<<"\n\n#Allele Frequencies Associated with this Penetrance Table\n";
	
	vector<string> genoLabels;
	BuildGenotypeLabels(genoLabels, loci.size());

	
	char label = 'A';
	char diff = 'a' - 'A';
	for (uint i=0; i<diseaseLoci.size(); i++) {
		Locus *l = diseaseLoci[i];
		file<<"FREQ "<<label<<" "<<l->Freq1()<<"\n";
		file<<"FREQ "<<label+diff<<" "<<l->Freq2()<<"\n";
		label++;
	}
	file<<"\n\n# The Penetrances are listed below\n";
	for (uint i=0; i<penCount; i++) 
		file<<genoLabels[i]<<"\t"<<penList[i]<<"\n";

}

void PenetranceModel::ReportDiseaseLoci(ostream &os, vector<Locus*>&diseaseLoci) {
	os<<"<TABLE border='1'><TR bgcolor=\"#dddddd\"><TH><Locus></TH><TH>Label</TH><TH>Freq Al. 1</TH><TH>Freq Al. 2</TH></TR>\n";
	
	char locus='A';
	for (uint i=0; i<diseaseLoci.size(); i++) {
		Locus *l = diseaseLoci[i];
		os<<"<TR><TD>"<<locus++<<"</TD>"
			<<"<TD>"<<l->GetLabel()<<"</TD>"
			<<"<TD>"<<l->Freq1()<<"</TD>"
			<<"<TD>"<<l->Freq2()<<"</TD></TR>";
	}
	os<<"</TABLE>";
}

void PenetranceModel::ReportPenetranceTable(ostream &os ){
	os<<"<P><CENTER><H4>Penetrance Table</H4><P>";
	os<<"<TABLE border='1'><TR bgcolor=\"#dddddd\"><TH>Cell Label</TH><TH>Penetrace</TH></TR>\n";

	vector<string> genoLabels;
	BuildGenotypeLabels(genoLabels, loci.size());

	for (uint i=0; i<penCount; i++) 
		os<<"<TR><TD>"<<genoLabels[i]<<"</TD><TD>"<<penList[i]<<"</TD></TR>\n";
	
	os<<"</TABLE></CENTER>\n";
	

}

uint PenetranceModel::GetGenotypeIdx(const char *genotype, uint lociCount) {
	// have to split each genotype every 2 letters is one genotype
	// A = 0 and a = 1 
	char currentConversion = 'A';
	uint index = 0;
	string gt(genotype);
	int pos = gt.find_first_of("Aa");


	for(uint currLoc=0; currLoc < lociCount; currLoc++)	{
		int gtValue = 0;

		char letter = genotype[pos++];
		if(letter != currentConversion) {
			gtValue = 1;
		}
		letter = genotype[pos++];
		if(letter == currentConversion)	{
			gtValue += 0;
		}
		else{
			gtValue += 1;
		}
		int multiplier = GetMultiplier(gtValue, currLoc);
		//cout<<genotype[pos-2]<<genotype[pos-1]<<" "<<multiplier<<" ";
		index += multiplier;
		
		currentConversion++;
	}
	return index;
}


string PenetranceModel::BuildGenotypeLabel(uint genotype, uint position) {
	char A='A'+position;
	char a='a'+position;
	
	char *label = new char[3];
	if (genotype == 0)
		sprintf(label, "%c%c", A, A);
	else if (genotype == 1)
		sprintf(label, "%c%c", A, a);
	else if (genotype == 2)
		sprintf(label, "%c%c", a, a);
	string finalLabel = label;
	delete[] label;
	return finalLabel;
}

string PenetranceModel::BuildGenotypeLabel(uint *genotypes, uint modelSize) {
	stringstream ss;
	for (uint i=0; i<modelSize; i++) 
		ss<<BuildGenotypeLabel(genotypes[i], i);
	return ss.str();
}



void PenetranceModel::BuildGenotypeLabels(vector<string>& labels, uint modelSize) {
	uint *genotypes = new uint[modelSize+1];
	uint position = 0;

//	for (uint i=0; i<modelSize; i++) 
//		genotypes[i]=0;
	memset((void*)genotypes, 0, (modelSize + 1 )*sizeof(uint));

	labels.clear();	
	position = modelSize - 1;

	while (genotypes[0]<3) {	
		string label = BuildGenotypeLabel(genotypes, modelSize);
		
		labels.push_back(label);
		//cout<<id++<<"\t"<<label<<"\n";

		if (++genotypes[position]>2 && position > 0) {
			//Find the highest position of rollover
			while (position-- > 0 && ++genotypes[position] > 2) {}

			while (position < modelSize - 1) 
				genotypes[++position] = 0;
		}
	}	
	delete[] genotypes;
}

/**	
 * @brief We are assuming that idx will be 0 based and ordered correctly
 * AABBCCDD, AABBCCDd, AABBCCdd, AABBCcDD, AABBCcDd, etc....
 */
void PenetranceModel::AddPenetrance(uint idx, double penetrance) {
	assert(idx < penCount && (int)idx > -1);
	penList[idx]=penetrance;
}

PenetranceModel::~PenetranceModel() {
	if (totalQueries > 0)
		cout<<"Disease model - Observed Prevalence: "<<(double)observed/(double)totalQueries<<"\n";
	if (penList)
		delete[] penList;
}
DiseaseLocus PenetranceModel::GetLocus(istream& ss) {
	string word;							///<Temporary holder for a give word on the line
	DiseaseLocus l;
	ss>>word;
	int chrID = atoi(word.c_str());

	if (word.length() == 0)
		return l;
#ifdef USE_XY
	if (chrID == 0) {
#else
	if (chrID < 1) {
#endif //USE_XY
		string errMsg;
		errMsg = "Invalid Chromosome ID, " + word + ". Chromosome indices start with 1";
		throw Utility::Exception::General(errMsg.c_str());
	}
		
	if (!ss.eof()) {
		ss>>word;
		int locID = atoi(word.c_str());
		if (locID < 1) {
			string errMsg;
			errMsg = "Invalid Locus ID, " + word + ". Locus IDs must be valid numbers greater than 0\n";
			throw Utility::Exception::General(errMsg.c_str());
		}
#ifdef USE_XY
		if (chrID == -1) {
			l.chromosome = -1;
		}
		else 
#endif //USE_XY
			l.chromosome = chrID - 1;
		l.locusIdx   = locID - 1;
	}
	return l;
}

int PenetranceModel::BuildGenotypeIndex(int *genotypes, int numGenos) {
    int index=0;
    for(int i=0; i<numGenos; i++){
		index += GetMultiplier(genotypes[i], i);
    }

	return index;	
}

int PenetranceModel::ConvertGenotype(int index, int gt) {
	return gt;
}

int PenetranceModel::BuildGenotypeIndex(vector<uint> & genotypes) {
    uint numGenos = genotypes.size();
    uint index=0;


    for(uint i=0; i<numGenos; i++){
#ifdef DEBUG_PRODUCTION
		cout<<genotypes[i]<<" ";
#endif
		index += GetMultiplier(genotypes[i], i);
    }
#ifdef DEBUG_PRODUCTION
		cout<<" == "<<index<<"\n";
#endif

	return index;	
}

int PenetranceModel::GetStatus(std::vector<uint>& genotypes, float& outcome) {
	totalQueries++;
#ifdef DEBUG_PRODUCTION
	cout<<"IsAffected: ";
#endif
	uint index = BuildGenotypeIndex(genotypes);

	assert(index <penCount);
	double risk = penList[index];
	double draw = Random::globalGenerator.drand();
	if (draw < risk) {
		observed++;
		outcome = 1.0;
		return 1;
	}
	outcome = 0.0;
	return 0;
}


PenetranceModel *PenetranceModel::GetModel(istream& ss, PoolManager *pools) {
	PenetranceModel *model = NULL;
	string name;
	string cmd;
	string modelType;					///<Simpen, Penetrance or simla

	ss>>modelType;			//>>percentage;

	if (strcmp(modelType.c_str(), PenTable) == 0) {
		model = new PenFileModel();
	}	
	else if (strcmp(modelType.c_str(), Simpen) == 0) 
		model = new GAModel();
	else if (strcmp(modelType.c_str(), Simla) == 0)
		model = new SimlaModel();
	else if (strcmp(modelType.c_str(), Continuous) == 0)
		model = new ContinuousModel();
	if (model && !model->Init(ss, pools)) {
		delete model;
		model = NULL;
	}
	return model;
}

string PenetranceModel::ParseFilename(istream &s) {
	string filename;
	s>>filename; 

	size_t lastPos = filename.length() - 1;
	if (filename[0] == '"' && filename[lastPos] !='"') {
		//Parsing filename with spaces
		bool doContinue = true;
		while (doContinue && !s.eof()) {
			string word;
			s>>word;
			filename+=" " + word;
			doContinue = word[word.length() -1]!='"';
		}
		if (doContinue) {
			stringstream ss;
			ss<<"Unterminated quote was encountered in the configuration of the disease model\n"<<filename<<"\n";
			throw Utility::Exception::General(ss.str().c_str());
		}
	}
	return StripQuotes(filename.c_str());
}
/*
bool PenetranceModel::InvertAllele(vector<double>& table, int modelSize, int locusToInvert) {
	vector<int> genotypes;
	
	vector<double> tableCopy = table;
	
	totalGtCount = pow(3, modelSize - 1);
	for (int i=0; i<totalGtCount; i++) {
		ConvertToGenotypes(genotypes, modelSize, i);
		if (genotypes[locusToInvert - 1] == 0) 
			genotypes[locusToInvert - 1] = 2;
		else if (genotypes[locusToInvert - 1] == 2)
			genotypes[locusToInvert - 1] = 0);
		int destGt;
	
		ConvertFromGenotypes(genotypes, modelSize, destGt);
		table[destGt] = tableCopy[i];
	}
}

void PenetranceModel::ConvertFromGenotypes(vector<int> genotypes, int genotypeCount, int mlGenotype) {
	int currOffset = 1;
	
	mlGenotype = 0;;
	for (int i=0; i<genotypeCount; i++) {
		mlGenotype+=(currOffset*genotypes[i]);
		currOffset*=3;
	}
}

void PenetranceModel::ConvertToGenotypes(vector<int> genotypes, int genotypeCount, int mlGenotype) {
	genotypes.clear();
	genotypes.resize(genotypeCount);
	int curOffset = pow(3, genotypeCount-1);
	int remainderGT = mlGenotype;
	
	for (int i=genotypeCount - 1; i>=0; i--) {
		int gt = remainderGT/currOffset;
		remainderGT = remainderGT%currOffset;
		genotypes[i]=gt;
		curOffset /= 3;
	}	
}*/


}

}


