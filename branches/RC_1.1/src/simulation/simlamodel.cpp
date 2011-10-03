//
// C++ Implementation: simlamodel
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "simlamodel.h"

namespace Simulation {

namespace StatusModel {
/*
PREVALENCE 				0.05

LOCUS_COUNT 			1

# Since this is for single locus models, we can only have 1 locus models
MAX_INTEREACTION_SIZE	1

# Model type is as follows:
# 0 			- recessive
# 0 < N < 1.0 	- mutiplicative 
# 1				- dominant
# LOCUS_MODEL [locus idx] [model type] [Odds Ratio]
DEFINE_LOCUS			1 0 2

# For multilocus interactions, give the loci and the odds ratio
#ADD_INTERACTION 		1x2 
*/


SimlaModel::SimlaModel() : isInitialized(false) { }

SimlaModel::~SimlaModel() { }



bool SimlaModel::ConfigureLocus(uint chrID, uint locID, float beta, float dxFreq, float type, bool dxAtAl1) {
	char key[1024];
	sprintf(key, "%d:%d", chrID + 1, locID+1);

//	cout<<"Adding beta values for locus: "<<key<<" ("<<beta<<") with the allele frequency of: "<<maf<<"\n";
	ModelDetails &b = betaValues[key];
	b = ModelDetails(betaValues.size(), key, beta,type, dxFreq, dxAtAl1);
	cout<<"Adding beta values for locus("<<betaValues.size()<<"): "<<key<<" ("<<beta<<") with the allele frequency of: "<<dxFreq<<"\n";
	betaKeys.push_back(key);
	return true;
}
void SimlaModel::Refresh(PoolManager *pools) {
	size_t lociCount = loci.size();
	
	for (size_t i=0; i<lociCount; i++) {
		DiseaseLocus &locus = loci[i];
		pools->GetAlleleFrequency(locus.chromosome, locus.locusIdx, locus.alFreq1, locus.alFreq2);
	}
}

bool SimlaModel::Init(istream& ss, PoolManager *pools) {
	string word;							///<Used during parsing of the command line
	string labelType;						///<index or label
	bool success=true;
	
	loci.clear();
	betaKeys.clear();
	betaValues.clear();
	ss>>labelType;

	configFilename = ParseFilename(ss);

	//If the definition is wrong, let's abort
	if (labelType != "INDEX" && labelType != "LABEL")
		return false;
	
	prevalence = 0.0;

	poolMgr = pools;

	chromCount = pools->GetPoolCount();

	//Get prevelance
	ss>>prevalence>>interactionCount;

	if (prevalence > 0.0 && interactionCount > 0 ) {
		bool minIsDisease;
		while (!ss.eof()) {		
			float dx = 0.0, al1 = 0.0;
			DiseaseLocus loc;
			string diseaseAllele;
			if (labelType == "INDEX") {
				loc = GetLocus(ss);
				ss>>diseaseAllele;

				pools->GetAlleleFrequency( loc.chromosome, loc.locusIdx, al1, dx);

				if (diseaseAllele == "MAJ") {
					minIsDisease = false;
					if (al1 > dx) 
						dx = al1;
				}
				else if (diseaseAllele == "MIN" ) {
					minIsDisease = true;
					if (al1 < dx)
						dx = al1;
				}
				else 
					return loci.size() > 0;
			}
			else {
				string label;
				ss>>label;
				ss>>diseaseAllele;
				Locus l;
				if (pools->GetLocusReference(label.c_str(), l) ) {
					loc.chromosome = l.GetChromID();
					loc.locusIdx = l.GetID();
					al1 = l.Freq1();
					if (diseaseAllele == "MAJ")  {
						minIsDisease = false;
						dx = 1.0 - l.GetMinAlleleFreq();
					}
					else if (diseaseAllele == "MIN") {
						minIsDisease = true;
						dx = l.GetMinAlleleFreq();
					}
					 else 						
						return false;
				}
				//al2 = l.locus->GetAlleleFreq(1);
			}
#ifdef USE_XY
			if (loc.chromosome > 0 || loc.chromosome == -1) {
#else
			if (loc.chromosome >= 0) {
#endif
				float beta = 0, modelType = 0;
				ss>>beta>>modelType;			

				bool dxOnAllele1 = minIsDisease && (dx == al1) || !minIsDisease && dx == al1;
						
				if (ConfigureLocus(loc.chromosome, loc.locusIdx, beta, dx, modelType, dxOnAllele1))
					AddDiseaseLoci(loc.label.c_str(), loc.chromosome, loc.locusIdx, dx);
				else {
					cout<<"Unable to set beta value for locus: "<<loc.chromosome<<":"<<loc.locusIdx<<" ("<<word<<") Unable to continue\n";
					abort();
				}
			}
		}
		SetModelSize(loci.size());
	}
	else {
		cout<<"A Prevalence of "<<word<<" is invalid. Prevalence must be greater than 0.0\n";
		abort();	
	}
	return success;
}

void SimlaModel::SetupInteractions() {
	string model;
	string val;
	float beta = 0;

	int loci[64];

	map<string, ModelDetails>::iterator cur = betaValues.begin();
	map<string, ModelDetails>::iterator end = betaValues.end();
	
	for (; cur!=end; cur++ )  {
		ModelDetails &m = cur->second;
		int index = m.idx;
		penetranceModel.add_interaction(1, &index, m.beta);

		cout<<"\t\tAdd Beta for Main Effect "<<index<<" ("<<m.locusID<<"): "<<m.beta<<"\n";
	}
	
	if (strcmp(configFilename.c_str(), "NO_INTERACTIONS") == 0)
		return;

	ifstream file(StripQuotes(configFilename.c_str()).c_str(), ios::in);
	if(!file){
		cerr << "Simla Configuration File: " << configFilename << " can't be opened.\n";
		exit(1);
	}

	while (!file.eof()) {
		model="";
		file>>model>>val;
	
		if (model.length() > 0) {
			beta = atof(val.c_str());
			char *start = (char *)model.c_str();		//Start of a given number in string
			char *lastChar = start + strlen(start);
			char *end = strchr(start, 'x');	//End position of the number in the string
			
			char snpID[64];
					
			bool doContinue = true;
			int idx = 0;
			//Build the model 1 step at a time
			while (doContinue) {
				uint len = lastChar - start;
				if (end)
					len=end-start;
				strncpy(snpID, start, len);
				snpID[len]='\0';
				
				if ((size_t)atoi(snpID) > this->loci.size()) {
					cerr<<"\n\nInvalid SIMLA interaction configuration: "<<model<<" contains one or more locus definition that falls outside the number of model loci. SIMLA interactions are to be identified according to their index within the model itself. i.e. 1x3x4 means the first, 3rd and 4th locus of a 4 SNP model (or larger). These SNPS are in the order specified on the simla model definition. ";
					exit(1);
				}
				loci[idx++]=atoi(snpID);
				
				doContinue = end != NULL;
				if (doContinue) {
					start = end + 1;
					end = strchr(start, 'x');
				}
			}
			
			cout<<"\t\tAdd Interaction: \t";
			for (int i = 0; i<idx; i++) 
				cout<<loci[i]<<" ";
			cout<<"\tBeta: "<<beta<<"\n";
	
			penetranceModel.add_interaction(idx, loci, beta);
		}
	}	
}


string SimlaModel::GetModelConfiguration() {
	static string minmax[]={"MIN","MAJ"};
	stringstream ss;
	ss<<"DEFINE_MODEL SIMLA INDEX \""<<configFilename<<"\" "<<prevalence<<" "<<interactionCount;
	map<string, ModelDetails>::iterator itr = betaValues.begin();
	map<string, ModelDetails>::iterator end = betaValues.end();

	while (itr != end) {
		ModelDetails &details=itr->second;
		cout<<"Model Details: "<<details.locusID<<" "<<details.dxAtAllele1<<" "<<details.beta<<" "<<details.modelType<<"\n";

		ss<<"    ";

		size_t posCol = details.locusID.find(':', 0); 
		if (posCol != string::npos)
			ss<<details.locusID.replace(posCol, 1, " ");
		else
			ss<<details.locusID;
		
		ss<<" "<<minmax[details.dxAtAllele1]<<" "<<details.beta<<" "<<details.modelType;
		itr++;
	}
	return ss.str();
}

void SimlaModel::Load() {
	if (isInitialized)  {
		cout<<"Cowardly refusing to reinitialize Simla with the configuration: "<<configFilename<<"\n";
		return;
	}

	uint modelSize = loci.size();
	penetranceModel.start_model(modelSize, interactionCount, prevalence);
	cout<<"Simla Model Initialization Begun: \n\tLoci Count: "<<modelSize<<"\n\tMax Interaction Size: "<<interactionCount<<"\n\tPrevalence: "<<prevalence<<"\n";


	assert(betaKeys.size() == modelSize);
	for (uint i=0; i<modelSize; i++) {
		ModelDetails &m = betaValues[betaKeys[i]];
		penetranceModel.add_locus(m.modelType, m.dxFreq);
		cout<<"\t\tAdd Locus: "<<m.locusID<<"\t"<<m.modelType<<"\t"<<m.dxFreq<<"\n";
	}
	
	SetupInteractions();

	uint seedSize = penetranceModel.close_model();

	int *seedData = new int[seedSize];
	
	GenerateSeedData(seedSize, seedData);
	//penetranceModel.find_x0(seedData);
	isInitialized = true;


	cout<<"\nDisease model in use: "<<configFilename<<"\n";
	int headerWidth = 45;
	//cout<<setw(4)<<setprecision(3)<<(100.0 * probability)<<"% \n";
	
	cout<<setw(headerWidth-15)<<""<<loci[0].label;
	for (uint i=1; i<modelSize; i++) 
		cout<<"x"<<loci[0].label;
	cout<<endl;
		
	cout<<"          --Penetrance table: \n";
	char geno = 'A';
	for (uint i=0; i<loci.size(); i++) 
		cout<<"\tFreq("<<(char)(geno+i)<<"): "<<loci[i].alFreq1<<"\n";

	vector<string> pentable;
	BuildGenotypeLabels(pentable, modelSize);
	for (uint i=0; i<penCount; i++)
		cout<<"              "<<pentable[i]<<" : "<<penList[i]<<"\n";

	delete[] seedData;
}


void SimlaModel::GenerateReport(ostream& os, uint headerWidth) {
	os<<setw(headerWidth)<<right<<"Model Type: "<<Details()<<"\n";
}

void SimlaModel::GenerateDetailedReport(ostream &os, vector<Locus*> &diseaseLoci) {
	os<<"<H3>SIMLA Based Model:</H3>\n";
	if (loci.size() != diseaseLoci.size()) 
		os<<"<B>Invalid Locus List Specified</B>";

	ReportPenetranceTable(os);
	
}

int SimlaModel::BuildConvertedGenotypeIndex(vector<uint> & genotypes) {
    uint numGenos = genotypes.size();
    uint index=0;


    for(uint i=0; i<numGenos; i++){
#ifdef DEBUG_PRODUCTION
		cout<<genotypes[i]<<" ";
#endif
		index += GetMultiplier( ConvertGenotype(i, genotypes[i]), i);
    }
#ifdef DEBUG_PRODUCTION
		cout<<" == "<<index<<"\n";
#endif

	return index;	
}


int SimlaModel::ConvertSimlaGenotype(int mlGenotype) {
	int modelSize = loci.size();

	int locus = 0, index = 0;
	int curGenotype = mlGenotype;
	for (int i = modelSize - 1; i> 0; i--) {
		int div = (int)powf(3.0, i);
		int localGt = curGenotype / div;
		curGenotype = curGenotype % div;
		index+=(div * ConvertGenotype(locus++, localGt));
	}
	index+=ConvertGenotype(locus++, curGenotype);
	return index;
}

void SimlaModel::GenerateSeedData(int seedSize, int *seedData) {
	size_t modelSize = loci.size();
	for (int i=0; i<seedSize; i++) {
		//Eventually we will avoid messing with the whole genome
		Individual ind(0, 0, chromCount);
		DrawIndividual(ind);

		vector<uint> genotypes;
		for (uint p=0; p<modelSize; p++)  
			genotypes.push_back(ind.GetGenotype(loci[p].chromosome, loci[p].locusIdx));

#ifdef DEBUG_PRODUCTION
		cout<<"GenerateSeedData ("<<i<<")\t";
#endif	
		int seedValue = BuildConvertedGenotypeIndex(genotypes);
		seedData[i] = seedValue;

	}
	penetranceModel.find_x0(seedData);

	penCount = (uint)powf(3.0, (modelSize));
	penList = new double[penCount];
	
	//OK, for now, we'll capture each of the possible genotypes. 
	//This might be better to do one at a time, but this seems likely to be better
	for (uint i=0; i<penCount; i++) {
		penList[ConvertSimlaGenotype(i)] = penetranceModel.fx(i);

#ifdef DEBUG_PRODUCTION
		int l[modelSize]; 
		cout<<"\t"<<i<<": "<<penList[ConvertSimlaGenotype(i)]<<"\t(";
		penetranceModel.number2genotypes(i, l);
		for (uint m=0; m<modelSize; m++) 
			cout<<" "<<l[m];
		cout<<")\t"<<penetranceModel.get_index(l)<<"\t"<<BuildGenotypeIndex(l, modelSize)<<"\n";
#endif
	}

}

int SimlaModel::ConvertGenotype(int index, int genotype) {
	//For now, 2 is homozygous for minor allele with simla
	ModelDetails &details = betaValues[betaKeys[index]];
	static int conversions[] = { 0, 1, 2, 2, 1, 0 };
	
	if (!details.dxAtAllele1) 
		return conversions[genotype];
	
	else {
		return conversions[genotype+3];
	}
}


void SimlaModel::DrawIndividual(Individual &person) {
	poolMgr->DrawIndividual(person);

}

void SimlaModel::ReserveIndividual(Individual &person) {
	PoolManager::Iterator itr = poolMgr->GetIterator();
	ChromPool *pool = itr.GetNext();
	while (pool) {
		pool->ReserveIndividual( &person);
		pool  = itr.GetNext();
	}

}	

}

}
