//
// C++ Implementation: snpaligned
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "snpaligned.h"
#include "snppool.h"
#include <boost/tokenizer.hpp>

namespace Genetics {

using namespace boost;

SnpAligned::~SnpAligned()
{
	if (binaryLabel)
		delete[] binaryLabel;
	if (bitsets)
		delete[] bitsets;
	if (bitsetLabels)
		delete[] bitsetLabels;
}

SnpAligned::SnpAligned(const SnpAligned& other) {
	cout<<"An internal error has occured and running can't continue.\n";
	abort();
}


bool SnpAligned::SetValue(int indivId, uint value) {
//	asciiGenotypes[indivId]='0'+value;

	if (value > genotypeCount)	{
		/**
		 * @todo Work up a way to increase the number of genotypes
		 */
		cout<<"Too many genotypes attempted at this time!\n";
		assert(0);
	}

	BitSetType &gt=GetGenotype(value);
	assert((uint)indivId < gt.size());
	gt[indivId]=value > 0;
	return true;
}

/**
 * This simply cross &s the two genotype bitsets returning a new one that equates to the punett style matrix 
 */
SnpAligned *SnpAligned::punett(SnpAligned* other) {
	if (other == NULL)
		return NULL;
	SnpPool *pool=SnpPool::Instance();
	
	//When we cross, we don't cross the missing data (we merge them)
	uint genotypes=1 + (genotypeCount - 1) * (other->genotypeCount - 1);
	SnpAligned *combo = pool->GetSnp(genotypes);
	combo->Init(genotypes);
	uint otherGenos = other->genotypeCount;
	uint idx = 0;
	
	Label *comboLables = new Label[labelCount + other->labelCount];

	//Since Missing data shouldn't be considered at this time, Let's make sure that genotype(0) is propogated completely
	BitSetType &missingGenos=combo->GetGenotype( idx++ );	
	missingGenos=GetGenotype(0)&other->GetGenotype(0);
	string genotypeLabel;
	genotypeLabel=GetGenotypeLabel(0) + other->GetGenotypeLabel(0);
	combo->SetGenotypeLabel(0, genotypeLabel);
	string spacer = "   ";

	for (idx=0; idx<labelCount; idx++) 
		comboLables[idx]=binaryLabel[idx];
	for (uint i=0; i<other->labelCount; i++)
		comboLables[idx++]=other->binaryLabel[i];
	combo->SetLabel(idx, comboLables);

	idx = 1;
	//cout<<"Missing data: "<<missingGenos<<" <- "<<GetGenotype(0)<<" x "<<other->GetGenotype(0)<<"\n";
	for (uint i=1; i<genotypeCount; i++) {
		for (uint j=1; j<otherGenos; j++) {
			BitSetType &genos=combo->GetGenotype(idx++);
			genos=GetGenotype(i)&other->GetGenotype(j);
			genotypeLabel=GetGenotypeLabel(i) + spacer + other->GetGenotypeLabel(j);
			//cout<<"Punett "<<combo->GetLabel()<<" "<<genotypeLabel<<" - "<<idx - 1<<" "<<genos<<" <- "<<GetGenotype(i)<<" x "<<other->GetGenotype(j)<<endl;
			if (!genos.any()) {
				idx--;
				//combo->genotypeCount--;
			}
			else
				combo->SetGenotypeLabel(idx-1, genotypeLabel);
		}
	}
	combo->genotypeCount = idx;

	delete[] comboLables;
		
	pool->Release();
	return combo;
}

void SnpAligned::Reset() {
	ResetEvaluations();
	bitsets[0].set();
	for (uint i=1; i<genotypeCount; i++) 
		bitsets[i].reset();
}

string SnpAligned::toString() {
	stringstream strBuff;
	strBuff<<GetLabel()<<" ";
	for (uint i=0; i<individualCount; i++) {
		//Until we go to environmental data, we will use the boring approach
		if (bitsets[1][i])
			strBuff<<genLookup.GetHomozygote1()<<" ";
		else if (bitsets[2][i])
			strBuff<<genLookup.GetHeterozygote()<<" ";
		else if (bitsets[3][i])
			strBuff<<genLookup.GetHomozygote2()<<" ";
		else
			strBuff<<"MM";
	}
	return strBuff.str();
}

bool SnpAligned::Merge(const SnpAligned *other) {
	boost::dynamic_bitset<>::block_type itr;
	bool success=genotypeCount==other->genotypeCount;

	if (success) {
		individualCount+=other->individualCount;
		
		for (uint i=0; i<genotypeCount; i++) {
			to_block_range(other->bitsets[i],&itr);
			
			for (uint j=0; j<other->bitsets[i].num_blocks(); j++) 	{
				//cout<<bitsets[i]<<" + "<<other->bitsets[i]<<" = ";
				bitsets[i].append(itr++);
				//cout<<bitsets[i]<<"\n";	
			}
			bitsets[i].resize(individualCount,i>0);
		}
	}
	return success;
}

int SnpAligned::ImportBpSnp(uint linenumber, uint label, const char *line) {
	uint genoID = 0;
	uint genosFound = 0;
	int individuals=(1+strlen(line))/3;
	
	//Clear and resize the bitset
	if (individuals > (int)individualCount)
		InitBitsets(individuals);
		
	string tokenString = line;
	
	boost::char_separator<char> sep(" \t,\r", "", boost::drop_empty_tokens);;
	
	//typedef tokenizer<char_separator<char> > tokenizer;
	tokenizer<char_separator<char> > tok(tokenString, sep);
		
	tokenizer<char_separator<char> >::iterator i = tok.begin();
	tokenizer<char_separator<char> >::iterator end = tok.end();
	
	bitsets[0].set();
	//int pos=0;
	for (; i!=end; ++i) {
		genoID=genLookup.GetValue((*i).c_str());
		//asciiGenotypes[pos++]='0'+genoID;
		//Right now, using lance's approach, valid genotypes are 1-3
		if (genoID<4)
			GetGenotype(genoID)[genosFound++]=genoID > 0;
		else
			cout<<"Skipping genotype: "<<genoID<<".\n";
	}
//	if (pos>0)
//		asciiGenotypes[pos]='\0';

	SetLabel(Label(linenumber, label));
//	cout<<"Import completed("<<id<<"): "<<individuals<<" found: "<<asciiGenotypes<<"\n";
	return genosFound;
}

void SnpAligned::InitBitsets(uint individuals) {
	hrCellsSet = false;
	bitsets[0].resize(individuals, true);
	for (uint i=1; i < genotypeCount; i++){
		bitsets[i].resize(individuals);
		bitsets[i].reset();
	}

	individualCount=individuals;
}
int SnpAligned::ImportSnp(uint linenumber, uint label, const char *genotype) {
	int genCount = genotypeCount;
	uint individuals = strlen(genotype);
	int curGenotype;

	//Clear and resize the bitset
	InitBitsets(individuals);	
	bitsets[0].set();
	for (uint j=0; j<individuals; j++) {
		curGenotype = *(genotype+j) - '0';

		if (curGenotype > genCount)	{
			/**
			 * @todo Work up a way to increase the number of genotypes
	 		 */
			cout<<"Too many genotypes attempted at this time!\n";
			abort();
		}

		BitSetType &gt=GetGenotype(curGenotype);
		//Set the 0 bitsets to false since those should represent individuals that are present
		gt[j]=curGenotype>0;					
	}
	Label lbl(linenumber, label);
	SetLabel(lbl);
	return genCount;
}

int SnpAligned::ImportCompressedBinSnp(uint linenumber, uint label, uint values[], const int count/* = 1*/, const int framesize/* = 2*/)
{
	uint genotype=0;					///<This is used to determine the bit vector in which we will be storing the data


	/**
	 * @brief Help get the compressed binary figured out
	 */
	BinArrayParser<uint> binParser(values, count, framesize);
	bitsets[0].set();
	//int pos=0;
	for (uint i=0; i<individualCount; i++)		{
		genotype=binParser.GetNextValue();
		//asciiGenotypes[pos++]='0'+genotype;
		if (genotype<genotypeCount)	{
			BitSetType &gt=GetGenotype(genotype);
			gt[i]=genotype>0;
		}
	}
	//asciiGenotypes[pos]='\0';
	SetLabel(Label(linenumber, label));
	return individualCount;	
}

int SnpAligned::IsPerfect(const BitSetType &mask) {
	bool isPerfect = true;
	int perfectIndividuals=0;

	BitSetType t;
	for (uint i=1; isPerfect && i<genotypeCount; i++) {
		BitSetType aff   = bitsets[i] & mask;
		BitSetType unaff = bitsets[i] & ~mask;
	
		int affected=aff.count();
		int unaffected=unaff.count();
			
		if ((affected * unaffected)==0)
			perfectIndividuals+=abs(affected-unaffected);
		else
			isPerfect=false;
	}
	if (isPerfect)
		return perfectIndividuals;
	else	
		return 0;
}

int SnpAligned::IsPerfect(uint genotype, const BitSetType &mask) {
	BitSetType aff   = bitsets[genotype] & mask;
	BitSetType unaff = bitsets[genotype] & ~mask;

	int affected=aff.count();
	int unaffected=unaff.count();

	int perfectIndividuals=0;
	
	if ((affected * unaffected)==0)
		perfectIndividuals=affected-unaffected;
	return perfectIndividuals;		
}



}
