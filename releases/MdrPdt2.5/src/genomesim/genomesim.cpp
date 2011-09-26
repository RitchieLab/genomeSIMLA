//
// C++ Implementation: genomesim
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "genomesim.h"
#include <iostream>
#include <iomanip>
#include "utility/lineparser.h"


#ifndef TEST_APP
int main(int argc, char **argv) {
	GenomeSIM::GenomeSim app;

	if (app.ParseCmdLine(argc, argv))
		app.Start();
	return 0;
}
	
#endif

namespace GenomeSIM {

using namespace std;

GenomeSim::GenomeSim() : startGeneration(0), additionalGenerations(0) {
}


GenomeSim::~GenomeSim() {
}


/**
 * @brief Prints the help contents
 */
void GenomeSim::PrintHelp() {
	PrintBanner();
#ifdef USE_MPI
	cout<<"usage: genomeSIMp <configuration file> <project name> options\n";
#else
	cout<<"usage: genomeSIMs <project name> options\n";
#endif
	cout<<"\t-p (--project) name  : Specify part of the filenames generated\n";
	cout<<"\t-l (--load-generation) number : Indicates a specific generation to be loaded and run from\n";
	cout<<"\t-w (--write-datasets) on/off : Indicates that datasets should or should not be produced from the final execution\n";
	cout<<"\t-d (--drop-points) first freq count : Override the drop point configuration\n";
	cout<<"\nThis application is designed to produce simlulated genetic data\n";
}


/**	
 * @brief Starts execution
 */
void GenomeSim::Start() {

	enum Mode {
		InitialMode, LoadAndRun, LoadAndSave };

	//For now, let's just get the first part working
	PoolManager *pools=configuration.GetPoolManager();
	uint currGen = pools->GetCurrentGeneration();


	//Perform some basic reporting
	configuration.GenerateReport(cout);

	//OK, the first drop defined specifically	
	int moreGenerations =configuration.generalSettings.firstDropPoint - currGen;

	//Work out the drop schedule for the others
	int dropPoints = configuration.generalSettings.dropCount;
	uint dropFrequency = configuration.generalSettings.dropFrequency;

	cout<<"Chromosome Pools initialized. \n\t"<<setw(35)<<"Current Generation: "<<currGen<<"\n";
	if (moreGenerations > 0) {
		cout<<"Chromosome Pools initialized and ready to go. \n";
		cout<<"Running from generation "<<currGen<<" to generation "<<configuration.generalSettings.firstDropPoint<<"\n";
		pools->AdvanceGenerations(moreGenerations, configuration.GetGrowthRate());
		if (configuration.generalSettings.createLDMaps)
			pools->GenerateLDMaps(configuration.generalSettings.outputName.c_str(), configuration.generalSettings.keepPhasedOutput);
		else if (configuration.generalSettings.keepPhasedOutput)
			pools->Dump(configuration.generalSettings.outputName.c_str());
	}
	else if (dropFrequency > 0 && moreGenerations < 0) {
		//figure out how many drop points are left
		dropPoints += (int)(((float)moreGenerations)/(float)dropFrequency);
	}
	
	for (int i=0; i<dropPoints; i++) {
		pools->AdvanceGenerations(dropFrequency, configuration.GetGrowthRate());
		if (configuration.generalSettings.createLDMaps)
			pools->GenerateLDMaps(configuration.generalSettings.outputName.c_str(), configuration.generalSettings.keepPhasedOutput);
		else if (configuration.generalSettings.keepPhasedOutput)
			pools->Dump(configuration.generalSettings.outputName.c_str());
		cout<<pools->Size()<<" chromosomes written to file for generation: "<<pools->GetCurrentGeneration()<<"\n";
	}

	cout<<"Generation advancement completed. \n";
	cout<<"\tNumber of generational drops performed: "<<dropPoints<<"\n";
	cout<<"\tCurrent Generation:                     "<<pools->GetCurrentGeneration()<<"\n";
	

	if (configuration.datasetSettings.writeDatasets)  {
		cout<<"Generating datasets\n";
		//configuration.statusSettings.models.GenerateReport(cout, 45);
		cout<<"Writing "<<configuration.datasetSettings.samples.size()<<" samples (of "<<configuration.datasetSettings.simSets<<" files each)\n";
		WriteDatasets(configuration.datasetSettings.individualCount);
	}

}

void GenomeSim::WriteDatasets(uint indCount) {
	PoolManager *pools = configuration.GetPoolManager();
	ModelManager &models = configuration.statusSettings.models;

	models.LoadModelData();
	
	uint sampleCount = configuration.datasetSettings.samples.size();
	uint datasetCount = configuration.datasetSettings.simSets;
	for (uint i=0; i<sampleCount; i++) {
		cout<<"Writing datasets for sample "<<i<<"\n";
		for (uint p=0; p<datasetCount; p++) {
			cout<<"*";
			cout.flush();
			configuration.datasetSettings.samples[i]->BuildSample(*pools, models, indCount);
			configuration.datasetSettings.samples[i]->ApplyErrors(pools, indCount);


			int snpCount = configuration.datasetSettings.samples[i]->GetSnpCount();
			uint gtCount[snpCount*4];
			for (int n=0; n<snpCount; n++) {
				gtCount[n*4+0]=0;
				gtCount[n*4+1]=0;
				gtCount[n*4+2]=0;
				gtCount[n*4+3]=0;
			}



			stringstream filename;
			filename<<configuration.generalSettings.outputName;
			filename<<"."<<pools->GetCurrentGeneration();
			filename<<"."<<i+1<<"."<<p+1;
			filename<<configuration.datasetSettings.samples[i]->GetDescriptor();
			ofstream file(filename.str().c_str(), ios_base::out);
			int observations = configuration.datasetSettings.samples[i]->WriteDataset(file, gtCount);
			file.close();

			stringstream report;
			report<<filename.str()<<".report";
			file.open(report.str().c_str(), ios::out);
			file<<"\n\n";
			file<<"snp #, Total, ??, AA, Aa, aa\n";
			for (int n=0; n<snpCount; n++) {
				file<<i<<","<<observations<<","
					<<gtCount[n*4]<<","
					<<gtCount[n*4+1]<<","
					<<gtCount[n*4+2]<<","
					<<gtCount[n*4+3]<<endl;
			}
			file.close();

			if (configuration.generalSettings.createSampleLDmaps) {
				filename<<".hap";
				file.open(filename.str().c_str(), ios_base::out);
				configuration.datasetSettings.samples[i]->WritePhased(file);
				file.close();

				if (configuration.datasetSettings.samples[i]->GetSnpCount() > PoolManager::haploWindowSize)
					pools->GenerateLDMapOverview(filename.str().c_str(), NULL, true);
				else
					pools->GenerateLDMap(filename.str().c_str(), NULL, true);
			}
		}
		cout<<"\n";
	}
}


bool GenomeSim::GetBoolean(const char *value) {
	bool isOn=false;
	isOn = strcmp(value, "ON")==0 || strcmp(value, "on")==0 || strcmp(value, "On")==0 || strcmp(value, "YES")==0 || strcmp(value, "Yes") ==0 || strcmp(value, "Yes")==0;
	return isOn;
}


int GenomeSim::ParseCmd(int curr, int argc, char **argv) {
	int nextCmd = curr+1;
	if (strcmp(argv[curr], "-p") == 0 || strcmp(argv[curr], "--project") ==0) 
		configuration.generalSettings.outputName = argv[nextCmd++];
	else if (strcmp(argv[curr], "-l")==0 || strcmp(argv[curr], "--load-generation")==0)
		configuration.generalSettings.firstGeneration = atoi(argv[nextCmd++]);
	else if (strcmp(argv[curr], "-w")==0 || strcmp(argv[curr], "--write-datasets")==0) {
		if (argc<nextCmd+1) {
			cout<<"-w (--write-datasets) must be followed by 1 value: Yes/No (or On/Off)\n";
			return -1;
		}
		configuration.datasetSettings.writeDatasets = GetBoolean(argv[nextCmd++]);
	}
	else if (strcmp(argv[curr], "-d")==0 || strcmp(argv[curr], "--drop-points")==0) {
		if (argc<nextCmd+3) {
			cout<<"-d (--drop-points) must be followed by 3 values: first_drop drop_frequency drop_count\n";
			return -1;
		}
		configuration.generalSettings.firstDropPoint = atoi(argv[nextCmd++]);
		configuration.generalSettings.dropFrequency = atoi(argv[nextCmd++]);
		configuration.generalSettings.dropCount = atoi(argv[nextCmd++]);
	}
	else {
		cout<<"Unknown argument: "<<argv[curr]<<"\n";
		return -1;
	}

	return nextCmd;
}

/**
 * Grab the pieces from the command line
 */
bool GenomeSim::ParseCmdLine(int argc, char **argv) {
	if (argc < 2) {
		PrintHelp();
		return false;
	}

	configuration.generalSettings.outputName = argv[1];


	Utility::LineParser lp;
	if (lp.Parse(argv[1], &configuration) > 0) {
		//Work out any other cmd line arguments
		int i =2;
		for (; i<argc && i>1;) {
			i=ParseCmd(i, argc, argv);
		}
		
		bool success = configuration.Validate();
		if (i<0) 
			cmdError=true;
		else if (!success)
			cout<<configuration.errorMsg;
		else {
			if (!configuration.PostLoad()) {
				cout<<"Unable populate the pool. Are you sure that your pool data exists?\n";
				cmdError = true;
			}
		}
	}
	else
		cmdError=true;
	return !cmdError && argc > 1;
}


}
