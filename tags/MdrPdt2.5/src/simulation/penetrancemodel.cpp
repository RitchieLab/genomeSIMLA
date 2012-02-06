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
#include "individual.h"
#include "utility/random.h"
//#include "defaults.h"

namespace Simulation {

using namespace std;
using namespace Utility;


PenetranceModel::PenetranceModel(uint modelID, uint modelSize, float prob) :DiseaseModel(modelID, prob), modelSize(modelSize) {
	//We can assume that the number of entries in the table are based on a power of 3 due to the 2allele model
	penCount = (uint)powf(3.0, (modelSize));
	penList = new double(penCount);
}
PenetranceModel::PenetranceModel(uint modelID, float prob) : DiseaseModel(modelID, prob), penCount(0), penList(NULL), modelLoci(NULL), modelSize(modelSize) { }

uint PenetranceModel::GetMultiplier(uint genotype, uint position) {
	assert(genotype < 3);
	uint pwr=modelSize - position - 1;
	return (uint)((float)genotype * powf(3.0, (float)pwr));
}

void PenetranceModel::GenerateReport( ostream &os, uint headerWidth) {
	os<<setw(headerWidth)<<right<<""<<"#"<<modelID<<" ";
	os<<setw(4)<<setprecision(3)<<(100.0 * probability)<<"% ";
	
/*	os<<setw(headerWidth)<<""<<modelLoci[0].chromosomeID<<":"<<modelLoci[1].locusID;
	for (uint i=1; i<modelSize; i++) 
		os<<modelLoci[0].chromosomeID<<":"<<modelLoci[1].locusID;
	os<<endl;

	os<<setw(headerWidth)<<""<<penList[0];
	for (uint i=1; i<penCount; i++)
		os<<":"<<penList[i];
*/
	os<<endl;
}

uint PenetranceModel::GetModelSize() {
	return modelSize;
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


//For now, just grab the different numbers on a line
void PenetranceModel::SetupDiseaseLoci(const char *line) {
	vector<LocusType> loci;
	
	stringstream ss(line);
	string ignore;
	ss>>ignore;				//We really don't care about the first word

	int value;
	int chrID;
	while (!ss.eof()) {
		value=0;
		chrID=0;
		ss>>chrID;
		if (chrID > 0)  {
			 if (!ss.eof()) {
				ss>>value;
				loci.push_back(LocusType(chrID, value));
				//cout<<chrID<<":"<<value<<" ";
			}
			else
				cout<<"Disease Loci configuration appears to not have chromosome IDs. Please correct this "<<chrID<<"\n";
		}
	}
	modelSize=loci.size();
	penCount = (uint)pow((double)3.0, (double)modelSize);
	penList = new double[penCount];
	modelLoci = new LocusType[modelSize];
		
	//I need a way to figure out how to figure out which chromosome a loci appears on
	for (uint i=0;i<modelSize; i++) {	
		modelLoci[i].chromosomeID = loci[i].chromosomeID - 1;
		modelLoci[i].locusID = loci[i].locusID - 1;
	}

	
	
}


bool PenetranceModel::ParseLine(const char *line, uint val) {
	bool validEntry = false;
	if (line[0] > 64 && line[0] < 123) {
		stringstream ss(line);
		string key;
		ss>>key;

		if (strcmp(key.c_str(), LABEL_DISEASELOCI) == 0)   {
			validEntry = true;
			SetupDiseaseLoci(line);
		}
		else if (strcmp(key.c_str(), LABEL_PENTABLE) == 0) {
			validEntry = true;
		}
		else {
			if (modelSize > 0) {
				float penetrance;
				ss>>penetrance;
				AddPenetrance(GetGenotypeIdx(key.c_str(), modelSize), penetrance);
				validEntry = true;
			}
		}
	}
	return validEntry;
}

void PenetranceModel::Load() {
	Load(filename.c_str());
}
	
void PenetranceModel::Load(const char *filename) {
	LineParser lp;
	this->filename = filename;
	uint validLines = lp.Parse(filename, this);
	
	if (validLines < 3 || validLines < 2 + penCount) {
		cout<<"Problems were encountered attempting to load the model file: "<<filename<<".\nPlease check that it is complete and ordered properly.\n";
		abort();
	}
	
	cout<<"\nDisease model in use: "<<filename<<"\n";
	int headerWidth = 45;
	cout<<setw(headerWidth)<<right<<""<<"#"<<modelID<<" ";
	cout<<setw(4)<<setprecision(3)<<(100.0 * probability)<<"% \n";
	
	cout<<setw(headerWidth-15)<<""<<modelLoci[0].chromosomeID+1<<":"<<modelLoci[0].locusID+1;
	for (uint i=1; i<modelSize; i++) 
		cout<<"x"<<modelLoci[i].chromosomeID+1<<":"<<modelLoci[i].locusID+1<<"\n";
	cout<<endl;
		
	cout<<"          --Penetrance table: \n";
	for (uint i=0; i<penCount; i++)
		cout<<"              "<<i<<" : "<<penList[i]<<"\n";

	cout<<endl;

}



/**
 * @brief Add a single locus and penetrance to the table
 * @param locus The snp locus being described
 * @param penetrance 
 * @return Returns the size of the model thus far

void PenetranceModel::AddDiseaseLoci(uint idx, uint chromID, uint locus) {
	assert(idx < modelSize);
	modelLoci[idx]=LocusType(chromID, locus);	
} */

/**	
 * @brief We are assuming that idx will be 0 based and ordered correctly
 * AABBCCDD, AABBCCDd, AABBCCdd, AABBCcDD, AABBCcDd, etc....
 */
void PenetranceModel::AddPenetrance(uint idx, double penetrance) {
	assert(idx < penCount);
	penList[idx]=penetrance;
}

PenetranceModel::~PenetranceModel() {
	if (penList)
		delete[] penList;

	if (modelLoci)
		delete[] modelLoci;
}

uint PenetranceModel::GetLocus(uint pos) {
	assert(pos<modelSize);
	return modelLoci[pos].locusID;
}

uint PenetranceModel::GetChromID(uint pos) {
	assert(pos<modelSize);
	return modelLoci[pos].chromosomeID;
}




bool PenetranceModel::IsAffected(std::vector<uint>& genotypes) {
    uint numGenos = genotypes.size();
	//assert((2 * modelSize) == numGenos);
	//string gt[]={"AA", "Aa", "aa"};
    uint index=0;

	//Otherwise, we will use the genotype to determine the risk
    for(uint i=0; i<numGenos; i++){
		//cout<<gt[genotypes[i]]<<" ";
		index += GetMultiplier( genotypes[i], i);
    }
	double risk = penList[index];
	double draw = Random::globalGenerator.drand();
	
	//cout<<"- "<<index<<"\tRisk: "<<risk<<"\tDraw: "<<draw;

/*	if (draw < risk)
		cout<<" - Affected\n";
	else 
		cout<<" - Unaffected\n";
	*/
	return draw < risk;
}


void PenetranceModel::SetFilename(const char *filename) {
	this->filename=filename;
}


}
