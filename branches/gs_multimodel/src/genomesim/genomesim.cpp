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
#include "utility/strings.h"
#include "utility/executionlog.h"
#include "timestamp.h"

#ifndef TEST_APP
int main(int argc, char **argv) {
	GenomeSIM::GenomeSim app;
	try {
		if (app.ParseCmdLine(argc, argv))
			app.Start();
	} catch (Utility::Exception::General const& e) {
		cout<<"Unhandled exception: \n";
		cout<<e.GetErrorMessage()<<"\n";
	}
	return 0;
}
	
#endif

namespace GenomeSIM {

using namespace std;

GenomeSim::GenomeSim() : startGeneration(0), additionalGenerations(0) {
	appname 	= "genomeSIMLA";
	appfunction = "Forward time simulation for SNP data";
	major	 	= APPMAJOR;
	minor 	 	= APPMINOR;
	bugFixes 	= APPBUGFIX;
	buildNumber = BUILD_NUMBER;
	buildType	= BUILD_TYPE;
	buildDate	= BUILD_DATE;

	cout<<"-----------------------------------genomeSIMLA-------------------------------\n";
	cout<<appname<<" is a forward-time population simulation application intended to\n";
	cout<<"provide users the ability to produce datasets rich with LD patterns similar to\n";
	cout<<"those seen in real data.\n\n";
	cout<<"This version of genomeSIMLA is considered under beta. We highly recommend that you \n";
	cout<<"check the genomeSIMLA website to make sure that you are running the latest version: \n";
	cout<<"\thttp://chgr.mc.vanderbilt.edu/genomeSIMLA\n\n";
#ifdef DEBUG
	cerr<<"-------------------------------------DEBUG-----------------------------------\n";
	cerr<<"(This has been compiled for debugging, and might run considerably slower than in release model)\n";
	cerr<<"-----------------------------------------------------------------------------\n";
#endif

}


GenomeSim::~GenomeSim() {
}


/**
 * @brief Prints the help contents
 */
void GenomeSim::PrintHelp() {
	PrintBanner();
#ifdef USE_MPI
	cout<<"usage: genomeSIMLAp <configuration file> [datasets] [ld] <project name> options\n";
#else
	cout<<"usage: genomeSIMLAs <configuration file> [datasets] [ld] <project name> options\n";
#endif
	cout<<"\tdatasets - Activate dataset generations\n";
	cout<<"\tld       - Activate complete ld plot generation\n";
	cout<<"\t-p (--project) name  : Specify part of the filenames generated\n";
	cout<<"\t-l (--load-generation) number : Indicates a specific generation to be loaded and run from\n";
	cout<<"\t-d (--drop-points) first freq count : Override the drop point configuration\n";
	cout<<"\t-s (--seed) number : Override the seed value in the configuration file\n";

	cout<<"\nWhen running genomeSIMLA in either dataset or ld mode, no generational advancement will\n";
	cout<<"be performed. Therefore, you should be sure to load the appropriate generation using the\n";
	cout<<"-l flag. The line below will load generation 100 and generate detailed LD plots from it.\n";
	cout<<"\tgenomeSIMLAs config.sim ld -l 100\n";

	cout<<"\nIf you have questions or issues regarding genomeSIMLA, please direct them to the \n";
	cout<<"genomeSIMLA distribution list, genomeSIMLA@chgr.mc.vanderbilt.edu\n";
}


string GenomeSim::OpenSummaryReport(ofstream &summary) {
	string reportFileName;
	char *filename = new char[4096];

	uint generation = configuration.GetPoolManager()->GetCurrentGeneration();
	string project = configuration.generalSettings.outputName.c_str();
	string sampled = "";
	if (ChromPool::fastLD)
		sampled = "-sampled";
	sprintf(filename, "%s.Index%s-%d.ld.html", project.c_str(), sampled.c_str(), generation);
	reportFileName = filename;
	summary.open(filename);
	summary<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
	summary<<"<HTML><HEADER><TITLE> GenomeSIMLA Summary Report - "<<project<<" Generation "<<generation<<"</TITLE></HEADER>\n";
	summary<<"<BODY>\n";
	summary<<"<link rel=\"stylesheet\" type=\"text/css\" href='"<<ChromPool::cssFilename<<"'>\n<DIV class='image-frame'>\n";
	summary<<"<CENTER>\n<H1>GenomeSIMLA Summary Report</H1>\n";
	summary<<"\t<TABLE class='simple'><CAPTION>Execution Details</CAPTION>\n";
	summary<<"\t\t<TR><TH>Configuration</TH><TH>Project Name</TH><TH>Pool Size</TH></TR>\n";
	summary<<"\t\t<TR><TD>"<<configuration.GetConfigurationFilename()<<"</TD><TD>"<<project<<"</TD><TD>"<<configuration.GetGrowthRate()->GetPopulationSize(generation)<<"</TD></TR>\n\t</TABLE>\n";
	summary<<"\t<TABLE><CAPTION>Chromosomes</CAPTION>\n";

	delete[] filename;
	return reportFileName;
}

/**	
 * @brief Starts execution
 */
void GenomeSim::Start() {
	ExecutionLog log(configuration.GetConfigurationFilename().c_str());
	//For now, let's just get the first part working
	PoolManager *pools=configuration.GetPoolManager();
	uint currGen = pools->GetCurrentGeneration();
	string growthChart = configuration.GetGrowthChartFilename();

	log.InitProject(configuration.generalSettings.seed, currGen, configuration.generalSettings.outputName.c_str());

	//Perform some basic reporting
	configuration.GenerateReport(cout);

	//OK, the first drop defined specifically	
	int moreGenerations =configuration.generalSettings.firstDropPoint - currGen;

	//Work out the drop schedule for the others
	int dropPoints = configuration.generalSettings.dropCount;
	uint dropFrequency = configuration.generalSettings.dropFrequency;

	//If we want to write datasets or detailed LD, we don't want to advance anymore
	if (configuration.generalSettings.writeDetailedLD || configuration.datasetSettings.writeDatasets)
		dropPoints = moreGenerations = 0;

	cout<<"Chromosome Pools initialized. \n\t"<<setw(35)<<"Current Generation: "<<currGen<<"\n";
	string mainReportFilename = "No_Report_Created";
	string locusReportFilename = "No_Report_Created";
	if (moreGenerations > 0) {
		cout<<"Chromosome Pools initialized and ready to go. \n";
		cout<<"Running from generation "<<currGen<<" to generation "<<configuration.generalSettings.firstDropPoint<<"\n";
		log.StartRun(currGen+moreGenerations);
		pools->AdvanceGenerations(moreGenerations, configuration.GetGrowthRate(), configuration.generalSettings.closePoolsBetweenAdvancing);
		if (PoolManager::generateOverviewLD) {
			cout<<"Writing Overview LD\n";
			ofstream summary;
			mainReportFilename = OpenSummaryReport(summary);
			locusReportFilename = pools->ProduceOverviewLD(configuration.generalSettings.outputName.c_str(), false, configuration.generalSettings.closePoolsBetweenAdvancing, summary, growthChart.c_str());
		}
		log.AppendReport("Locus_Search", locusReportFilename.c_str());
		log.EndRun(pools->GetExpressionCount(), mainReportFilename.c_str(), locusReportFilename.c_str());

/*		if (ChromPool::produceHaploviewLDInput) {
			cout<<"Writing Haploview enabled LD Input:\n";
			pools->WriteLdInputForHaploview(configuration.generalSettings.outputName.c_str(), configuration.generalSettings.closePoolsBetweenAdvancing);
		}
*/
	}
	else if (dropFrequency > 0 && moreGenerations < 0) {
		//figure out how many drop points are left
		dropPoints += (int)(((float)moreGenerations)/(float)dropFrequency);
	}

	string fastLD = "";
	if (ChromPool::fastLD)
		fastLD="(Sampled)";
	for (int i=1; i<dropPoints; i++) {
		log.StartRun(currGen+moreGenerations+(dropFrequency*i));
		pools->AdvanceGenerations(dropFrequency, configuration.GetGrowthRate(), configuration.generalSettings.closePoolsBetweenAdvancing);
	
		string summaryReportFilename = "No_Report_Found";	
		string locusReportFilename = "No_Report_found";
		if (PoolManager::generateOverviewLD) {
			cout<<"Writing Overview LD\n";
			ofstream summary;
			summaryReportFilename = OpenSummaryReport(summary);
			locusReportFilename = pools->ProduceOverviewLD(configuration.generalSettings.outputName.c_str(), false, configuration.generalSettings.closePoolsBetweenAdvancing, summary, growthChart.c_str());
		}
		log.AppendReport("Locus_Report_Filename", locusReportFilename.c_str());
		log.EndRun(pools->GetExpressionCount(), summaryReportFilename.c_str(), (string("Index")+fastLD).c_str());
/*		if (ChromPool::produceHaploviewLDInput) {
			cout<<"Writing Haploview enabled LD Input:\n";
			pools->WriteLdInputForHaploview(configuration.generalSettings.outputName.c_str(), configuration.generalSettings.closePoolsBetweenAdvancing);
		}
*/
		cout<<pools->Size()<<" chromosomes written to file for generation: "<<pools->GetCurrentGeneration()<<"\n";
	}

	cout<<"Generation advancement completed. \n";
	cout<<"\tNumber of generational drops performed: "<<dropPoints<<"\n";
	cout<<"\tCurrent Generation:                     "<<pools->GetCurrentGeneration()<<"\n";
	

	if (configuration.datasetSettings.writeDatasets)  {
		cout<<"Generating datasets\n";
		//configuration.statusSettings.models.GenerateReport(cout, 45);
		cout<<"Writing "<<configuration.datasetSettings.samples.size()<<" samples (of "<<configuration.datasetSettings.simSets<<" files each)\n";
		WriteDatasets();
	}

	if (configuration.generalSettings.writeDetailedLD) {
		ChromPool::fastLD = false;
		ofstream summary;
		OpenSummaryReport(summary);

		pools->ProduceOverviewLD(configuration.generalSettings.outputName.c_str(), false, configuration.generalSettings.closePoolsBetweenAdvancing, summary, growthChart.c_str());
	}


	cout<<"\nDone!\n";

}

void GenomeSim::WriteDatasets() {
	stringstream fileprefix;
	PoolManager *pools = configuration.GetPoolManager();


	size_t popSize = pools->GetExpressionCount();
	if (popSize < 8000) {
		cout<<"Error. Pool size is too small. Please adjust the growth rate so that the pool from which the datasets are generated has at least 8,000 individuals. We recommend populations of over 20K, when possible- and even higher in the event of extremely rare diseases.\n";
		abort();
	}
	else if (popSize < 20000) 
		cout<<"Warning!\nThe current population is: "<<popSize<<". We highly recommend large populations to ensure as little redundant genetic information as possible. \n";

	PenetranceModel *model = configuration.statusSettings.LoadModel(configuration.generalSettings.outputName.c_str(), pools);

	fileprefix<<configuration.generalSettings.outputName<<"."<<pools->GetCurrentGeneration();
	
	uint sampleCount = configuration.datasetSettings.samples.size();
	uint datasetCount = configuration.datasetSettings.simSets;
	
	string sampleMarkerInfo;

	cout<<"Model use: "<<model->Details()<<"\n";

	//Iterate through the various samples that are to be created
	for (uint i=0; i<sampleCount; i++) {
		if (!configuration.datasetSettings.samples[i]->Verify(model))
			abort();
		cout<<"\n\nWriting datasets for sample "<<configuration.datasetSettings.samples[i]->Details()<<"\n";
		for (uint p=0; p<datasetCount; p++) {
			cout<<"*";
			cout.flush();
			if (p==0) 
				configuration.datasetSettings.samples[i]->PrepSample(pools, model);
			configuration.datasetSettings.samples[i]->BuildSample(*pools, model);
			configuration.datasetSettings.samples[i]->ApplyPhenocopyError(pools, model);			
			pools->ResolveGenotypes(configuration.generalSettings.outputName.c_str(), pools->GetCurrentGeneration(), configuration.datasetSettings.samples[i], model);
			configuration.datasetSettings.samples[i]->ApplyErrors(pools);

			int snpCount = configuration.datasetSettings.samples[i]->GetSnpCount();
			uint *gtCount = new uint[snpCount*4];
			for (int n=0; n<snpCount; n++) {
				gtCount[n*4+0]=0;
				gtCount[n*4+1]=0;
				gtCount[n*4+2]=0;
				gtCount[n*4+3]=0;
			}
			
			int observations = 0;
			ofstream file;
			stringstream filename;
			filename<<fileprefix.str()<<"."<<i+1<<"."<<p+1;
			filename<<configuration.datasetSettings.samples[i]->GetDescriptor();
			if (configuration.datasetSettings.binaryDatasets) {
				stringstream meta(filename.str().c_str());
				meta<<".meta";
				filename<<".bin";
				file.open(filename.str().c_str(), ios_base::out);
				ofstream metadata;
				if (configuration.datasetSettings.samples[i]->UsesMeta())
					metadata.open(meta.str().c_str(), ios_base::out);
				observations = configuration.datasetSettings.samples[i]->WriteBinaryDataset(metadata, file);
				file.close();
				metadata.close();
			} else {
				if (p==0) 
					configuration.datasetSettings.samples[i]->WriteSnpDetails(filename.str().c_str(), pools);
				observations = configuration.datasetSettings.samples[i]->WriteDataset(filename.str().c_str(), gtCount);
			}
			configuration.datasetSettings.samples[i]->AppendGenotypeCountsToReport(model->GetDiseaseLoci());
			//configuration.datasetSettings.samples[i]->ReportGenotypeCounts(model->GetDiseaseLoci());
			configuration.datasetSettings.samples[i]->Purge();
			delete[] gtCount;
		}

		cout<<"\n";
		configuration.datasetSettings.samples[i]->ReportGenotypeCounts(model->GetDiseaseLoci(), cout);
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
		if (nextCmd < argc) {
			configuration.generalSettings.outputName = argv[nextCmd++];
		}
		else {
			cout<<"-p (project name) must be followed by a string value\n";
			return -1;
		}
	else if (strcmp(argv[curr], "-l")==0 || strcmp(argv[curr], "--load-generation" )==0)
		if (nextCmd < argc)  {
			configuration.generalSettings.firstGeneration = atoi(argv[nextCmd++]);
			configuration.generalSettings.doLoad = true;
		}
		else {
			cout<<"-l (load generation) must be followed by an integer value\n";
			return -1;
		}
	else if (strcmp(argv[curr], "datasets") == 0) 
		configuration.datasetSettings.writeDatasets = true;
	else if (strcmp(argv[curr], "ld") == 0)
		configuration.generalSettings.writeDetailedLD = true;
	/*else if (strcmp(argv[curr], "-w")==0 || strcmp(argv[curr], "--write-datasets" )==0) {
		if (argc<nextCmd+1) {
			cout<<"-w (--write-datasets) must be followed by 1 value: Yes/No (or On/Off)\n";
			return -1;
		}
		configuration.datasetSettings.writeDatasets = GetBoolean(argv[nextCmd++]);
	}*/
	else if (strcmp(argv[curr], "-d")==0 || strcmp(argv[curr], "--drop-points")==0) {
		if (argc<nextCmd+3) {
			cout<<"-d (--drop-points) must be followed by 3 values: first_drop drop_frequency drop_count\n";
			return -1;
		}
		configuration.generalSettings.firstDropPoint = atoi(argv[nextCmd++]);
		configuration.generalSettings.dropFrequency = atoi(argv[nextCmd++]);
		configuration.generalSettings.dropCount = atoi(argv[nextCmd++]);
	}
	else if (strcmp(argv[curr], "-s")==0 || strcmp(argv[curr], "--seed")==0) {
		configuration.SetRandomSeed(atoi(argv[nextCmd++]));
	}		
	else if (strcmp(argv[curr], "-i")==0 || strcmp(argv[curr], "--init")==0) {
		string initType = argv[nextCmd++];
		if (initType == "ADAMEVE") 
			ChromPool::UseAdamEve = true;
		else if (initType == "EDEN")
			ChromPool::UseEden = true;
		
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


	Utility::LineParser lp('#');
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
