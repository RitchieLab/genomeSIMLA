//
// C++ Implementation: diseasemodeldetails
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "diseasemodeldetails.h"
#include <fstream>
#include "utility/exception.h"
#include "simulation/pedigreereferencesample.h"

namespace GenomeSIM {

namespace GUI {


DiseaseModelDetails::DiseaseModelDetails(const char *filename) : filename(filename) {}
DiseaseModelDetails::~DiseaseModelDetails()	{}


void DiseaseModelDetails::Load() {
	FileToMap data;
	data.Parse(filename.c_str());

	PostLoad(data);	
}

void DiseaseModelDetails::Save() {
	Save(filename.c_str());	
}

DiseaseModelDetails *DiseaseModelDetails::OpenDiseaseModel(const char *filename) {
	string extension = ExtractExtension(filename);
	string type="UNKNOWN";
	
	if (extension == "pen")
		type="PENTABLE";
	else if (extension == "simpen")
		type="SIMPEN";
	else if (extension == "simla")
		type="SIMLA";
	
	if (type != "UNKNOWN")
		return OpenDiseaseModel(type.c_str(), filename);
	return NULL;
}

DiseaseModelDetails *DiseaseModelDetails::OpenDiseaseModel(
		const char *type, const char *filename) {
	DiseaseModelDetails *model = NULL;

	string modelType = type;

	if (modelType == "PENTABLE") 
		model = new PenetranceModelDetails(filename);
	else if (modelType == "SIMPEN")
		model = new SimpenModelDetails(filename);
	else if (modelType == "SIMLA")
		model = new SimlaModelDetails(filename);
	return model;
}

PenetranceModelDetails::PenetranceModelDetails(const char *filename) : 
		DiseaseModelDetails(filename), freqThreshold(0.0) { }
PenetranceModelDetails::~PenetranceModelDetails() { }

void PenetranceModelDetails::LoadLoci(vector<string> &freqs) {
	map<string, double> freqCache;
	size_t allCount = freqs.size();
	locusCount = allCount / 4;

	for (size_t i=0; i<allCount; i+=2) 
		freqCache[freqs[i]]=atof(freqs[i+1].c_str());

	char al1[]="A";
	char al2[]="a";


	for (int i=0; i<locusCount; i++) {
		double f1=freqCache[al1];
		double f2=freqCache[al2];
		Freqs f(f1, f2, al1[0]);
		frequencies.push_back(f);
		al1[0]++;
		al2[0]++;
	}
}

void PenetranceModelDetails::PostLoad(FileToMap& data) {
	freqThreshold = data.GetDouble("FREQ_THRESHOLD");
	vector<string> freqs;
	if (data.GetLines("FREQ", freqs))
		LoadLoci(freqs);
	else
		locusCount = 0;

	LoadPenetrances(data);
}



string PenetranceModelDetails::BuildGenotypeLabel(uint genotype, uint position) {
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

string PenetranceModelDetails::BuildGenotypeLabel(uint *genotypes) {
	stringstream ss;
	
	for (int i=0; i<locusCount; i++)  {
		ss<<BuildGenotypeLabel(genotypes[i], i);	
		cout<<genotypes[i]<<" ";
	}
	cout<<"= "<<ss.str()<<"\n";
	return ss.str();
}

void PenetranceModelDetails::BuildGenotypeLabels() {
	uint *genotypes = new uint[locusCount+1];
	uint position = 0;

	memset((void*)genotypes, 0, (locusCount + 1 )*sizeof(uint));

	position = locusCount- 1;

	cout<<"Locus Count: "<<locusCount<<"\n";
	while (genotypes[0]<3) {	
		string label = BuildGenotypeLabel(genotypes);
		
		penetrances.push_back(Penetrance(label.c_str()));
		
		if (++genotypes[position]>2 && position > 0) {
			//Find the highest position of rollover
			while (position-- > 0 && ++genotypes[position] > 2) {}

			while (position < locusCount - 1) 
				genotypes[++position] = 0;
		}
	}	
}

void PenetranceModelDetails::LoadPenetrances(FileToMap& data) {
	penetrances.clear();
	BuildGenotypeLabels();

 	vector<Penetrance>::iterator itr = penetrances.begin();
	vector<Penetrance>::iterator end = penetrances.end();
	
	while (itr != end) {
		itr->pen = atof(data.GetLine(itr->label.c_str()).c_str());	
		itr++;
	}
}

string PenetranceModelDetails::GetConfigurationDetails(vector<Locus*> &loci) {
	stringstream ss;
	ss<<"DEFINE_MODEL PENTABLE LABEL \""<<filename<<"\" ";
	int locCount = loci.size();
	for (int i=0; i<locCount; i++) 
		ss<<loci[i]->GetLabel()<<" ";

	return ss.str();
}
void PenetranceModelDetails::Save(const char *filename) {
	ofstream file(filename);

	//Need to work out a comment section	
	file<<"# The Threshold will help protect genomeSIMLA from using a penetrance\n";
	file<<"# table with data whose allele frequencies vary too greatly from the population\n";
	file<<"# for which the table was designed.\n";
	file<<"FREQ_THRESHOLD "<<freqThreshold<<"\n";
	file<<"\n\n#Allele Frequencies Associated with this Penetrance Table\n";

	{
		vector<Freqs>::iterator itr = frequencies.begin();
		vector<Freqs>::iterator end = frequencies.end();
	
		char diff = 'a' - 'A';
		while (itr != end) {
			Freqs& freq = *itr;
			file<<"FREQ "<<freq.label<<" "<<freq.al1<<"\n";
			file<<"FREQ "<<freq.label+diff<<" "<<freq.al2<<"\n";
			itr++;		
		}
	}
	
	file<<"\n\n# The Penetrances are listed below\n";
	
  	vector<Penetrance>::iterator itr = penetrances.begin();
	vector<Penetrance>::iterator end = penetrances.end();
	while (itr != end) 
		file<<itr->label<<" "<<itr->pen<<"\n";
}



SimpenModelDetails::SimpenModelDetails(const char *filename) : DiseaseModelDetails(filename) { }
SimpenModelDetails::~SimpenModelDetails() { }

void SimpenModelDetails::PostLoad(FileToMap& data) {
	heritWeight 	= atof(data.GetLine("HERITWEIGHT").c_str());
	heritTarget 	= atof(data.GetLine("HERIT").c_str());
	
	varianceWeight 	= atof(data.GetLine("MARGWEIGHT").c_str());
	varianceTarget	= atof(data.GetLine("MARGVAR").c_str());
	
	orWeight 		= atof(data.GetLine("ODDSWEIGHT").c_str());
	orTarget		= atof(data.GetLine("ODDSRATIO").c_str());
	
	penTarget		= atof(data.GetLine("PENTARGET").c_str());
	gaSettings		= data.GetLine("GA_SETTINGS").c_str();
}


void SimpenModelDetails::WriteToFile(ostream& file, const char *key, vector<string>& values) {
	file<<key;
	
	int count=values.size();
	for (int i=0;i<count; i++) 
		file<<" "<<values[i];
	file<<"\n";
}

string SimpenModelDetails::WriteSimpenConfigFile() {
	stringstream ss;
	ss<<filename<<".simpen-cfg";
	std::ofstream file(ss.str().c_str());
	
	file<<"#Heritability\n";
	file<<"HERITWEIGHT "<<heritWeight<<"\n";
	file<<"HERIT       "<<heritTarget<<"\n";
	file<<"\n#Marginal Variance\n";
	file<<"MARGWEIGHT "<<varianceWeight<<"\n";
	file<<"MARGVAR    "<<varianceTarget<<"\n";
	file<<"\n#Odds Ratio\n";
	file<<"PENTARGET  "<<penTarget<<"\n";
	file<<"\n\n--------------The rest of the settings below is associated with the GA\n";

	//Now, we grab the rest of the stuff from the ga settings file
	map<string, int> localValues;
	localValues["HERITWEIGHT"] 	= 0;
	localValues["HERIT"] 		= 0;
	localValues["MARGWEIGHT"]	= 0;
	localValues["MARGVAR"]		= 0;
	localValues["ODDSWEIGHT"] 	= 0;
	localValues["ODDSRATIO"]	= 0;
	localValues["PENTARGET"]	= 0;

	FileToMap data;
	try {
		data.Parse(gaSettings.c_str());
	}
	catch (Utility::Exception::FileNotFound& e) {
		stringstream ss;
		ss<<"Base GA settings file, "<<gaSettings.c_str()<<", we unable to be opened- and the model was not able to be properly saved for execution. Please select a valid GA settings file";
		throw Utility::Exception::General(ss.str().c_str());
	}
	vector<string> contents;
	vector<string> keys=data.GetKeys();
	
	vector<string>::iterator itr=keys.begin();
	vector<string>::iterator end=keys.end();

	map<string, int>::iterator notFound = localValues.end();
	while (itr != end) {
		const char *key = (*itr).c_str();
		if (localValues.find(key) != notFound) {
			data.GetLines(key, contents);
			WriteToFile(file, key, contents);
		}
		itr++;
	}
	file.close();
	
	return ss.str();
}



string SimpenModelDetails::GetConfigurationDetails(vector<Locus*> &modelLoci) {
	string configFile = WriteSimpenConfigFile();
	stringstream ss;
	ss<<"DEFINE_MODEL SIMPEN LABEL \""<<configFile<<"\"";

	int locCount = modelLoci.size();
	for (int i=0; i<locCount; i++) 
		ss<<" "<<modelLoci[i]->GetLabel();
	return ss.str();
}

void SimpenModelDetails::Save(const char *filename) {
	ofstream file(filename);
	
	file<<"#Heritability\n";
	file<<"HERITWEIGHT "<<heritWeight<<"\n";
	file<<"HERIT       "<<heritTarget<<"\n";
	file<<"\n#Marginal Variance\n";
	file<<"MARGWEIGHT "<<varianceWeight<<"\n";
	file<<"MARGVAR    "<<varianceTarget<<"\n";
	file<<"\n#Odds Ratio\n";
	file<<"PENTARGET  "<<penTarget<<"\n";
	file<<"\n#GA Settings (filename containing the GA details\n";
	file<<"GA_SETTINGS "<<gaSettings<<"\n";
}


SimlaModelDetails::SimlaModelDetails(const char *filename) : DiseaseModelDetails(filename) { }
SimlaModelDetails::~SimlaModelDetails() { }

void SimlaModelDetails::PostLoad(FileToMap& data) {
	prevalence = atof(data.GetLine("SIMLA_PREVALENCE").c_str());
	vector<string> list;
	
	locusCount = atoi(data.GetLine("SIMLA_LOC_COUNT").c_str());
	
	if (data.GetLines("SIMLA_LOCUS", list)) {
		assert(locusCount == (int)(list.size() / 4));
		
		for (int i=0; i<locusCount; i++) {
			LocusOR locus(list[i*4].c_str()[0], atoi(list[i*4+1].c_str()) == 1, atof(list[i*4+2].c_str()), atof(list[i*4+3].c_str()));
			loci.push_back(locus);
		}
	}

	list.clear();
	if (data.GetLines("SIMLA_INTERACTION", list)) {
		int interactionCount = list.size()/2;
		for (int i=0; i<interactionCount; i++) {
			double oddsRatio = atof(list[i*2+1].c_str());
			Interaction interaction(list[i*2].c_str(), oddsRatio);
			interactions.push_back(interaction);
		}
	}
	
}	

string SimlaModelDetails::GetConfigurationDetails(vector<Locus*> &modelLoci) {
	stringstream ss;
	static string minMaj[] = {"MIN","MAJ"};

	//Setup the interaction file
	ss<<filename<<".interactions";
	string interactionFile = ss.str();
	ofstream intFile(interactionFile.c_str());
	int intCount = interactions.size();
	cout<<"Writing Simla Interaction file: \n";
	for (int i=0; i<intCount; i++)  {
		cout<<interactions[i].label<<" -> "<<interactions[i].ConvertToSIMLA()<<"\t"<<ToBeta(interactions[i].oddsRatio)<<"\n";
		intFile<<interactions[i].ConvertToSIMLA()<<" "<<ToBeta(interactions[i].oddsRatio)<<" "<<"\n";
	}
	intFile.close();

	ss.str("");
	ss<<"DEFINE_MODEL SIMLA LABEL \""<<interactionFile<<"\" "<<prevalence<<" "<<loci.size()<<"	";
	int locCount = modelLoci.size();
	for (int i=0; i<locCount; i++) 
		ss<<modelLoci[i]->GetLabel()<<" "
			<<minMaj[loci[i].diseaseAtMajor]<<" "
			<<ToBeta(loci[i].oddsRatio)<<" "
			<<loci[i].modelType<<" ";
	return ss.str();
}

void SimlaModelDetails::Save(const char *filename) {
	ofstream file(filename);

	file<<"#Disease Prevalence\n"<<prevalence<<"\n\n";

	file<<"#Markers assocaited with the disease\n";
	for (int i=0; i<locusCount; i++) {
		LocusOR &locus = loci[i];
		//We should convert this to beta values for SIMLA (we use Rel. Risk)
		cout<<"Writing Simla Details: "<<locus.label<<" "<<locus.diseaseAtMajor<<" "<<ToBeta(locus.oddsRatio)<<"\n";
		file<<"SIMLA_LOCUS "<<locus.label<<" "<<locus.diseaseAtMajor<<" "<<locus.oddsRatio<<"\n";
	}

	int intCount = interactions.size();
	file<<"\n\n#Interactions\n";
	for (int i=0; i<intCount; i++) {
		cout<<"Writing Interaction Details: "<<interactions[i].label<<" "<<interactions[i].oddsRatio<<"\n";
		file<<"SIMLA_INTERACTION "<<interactions[i].label<<" "<<interactions[i].oddsRatio<<"\n";	
	}
}

}

}
