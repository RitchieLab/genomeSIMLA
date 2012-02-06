//
// C++ Implementation: appcontroller
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "appcontroller.h"
#include <iostream>
#include <wx/msgdlg.h>
#include <wx/app.h>
#include <wx/filename.h>
#include "genomesimla.h"
#include "utility/exception.h"

#define SIMLOCK pthread_mutex_lock(&simLock);
#define SIMUNLOCK pthread_mutex_unlock(&simLock);

namespace GenomeSIM {

namespace GUI {

using namespace std;


StringArray AppController::gaSettings;

AppController::AppController() : execLog(NULL), pools(NULL), simIsRunning(false) {
	pthread_mutex_init(&simLock, NULL);

	InitUserResources();	

	/// Populate the gaSettings list
	FileToArray fta;
	LineParser lp;
	try {
		lp.Parse("gaSettings.conf", &fta);
	} catch (Utility::Exception::FileNotFound& e) {	}

	gaSettings = fta.strings;
	if (gaSettings.size() < 1) {
		wxAppTraits *traits = wxGetApp().GetTraits();
		wxFileName datadir;
		datadir.Assign(traits->GetStandardPaths().GetUserDataDir(), "tight.ga");
		gaSettings.push_back(datadir.GetFullPath().c_str());
	}
}

void AppController::InitUserResources() {
	wxGetApp().SetAppName("wxGenomeSIMLA");
	wxAppTraits *traits = wxGetApp().GetTraits();
	string userDir = traits->GetStandardPaths().GetUserDataDir().c_str();
	string appDir = traits->GetStandardPaths().GetResourcesDir().c_str();
	string userConfigDir = traits->GetStandardPaths().GetUserConfigDir().c_str();
	
	cout<<"User Dir: "<<userDir<<"\n"
		<<"App Dir : "<<appDir<<"\n"
		<<"Usr Cfg : "<<userConfigDir<<"\n";

	wxFileName gaCfgFile;
	gaCfgFile.Assign(_(userDir.c_str()), "gaSettings.conf");

	char cmd[9600];

	if (!gaCfgFile.DirExists()) { 
		cout<<"Creating directory: "<<userDir<<"\n";
		gaCfgFile.Mkdir();
	}

	if (!gaCfgFile.FileExists()) {
		
		cout<<"Creating default ga settings: "<<appDir.c_str()<<"/*.ga -> "<<userDir.c_str()<<"\n";
		
	
		wxFileName gaMaster;
		sprintf(cmd, "cp %s/*.ga  %s", appDir.c_str(), userDir.c_str());
		int rv = system(cmd);
	}
	wxFileName simLog;
	simLog.Assign(_(userDir.c_str()), "genomeSIMLA.err");
	freopen(simLog.GetFullPath(), "w", stderr);
	cout<<"Error Log: "<<simLog.GetFullPath()<<"\n";
	simLog.Assign(_(userDir.c_str()), "genomeSIMLA.log");
	cout<<"Std Out: "<<simLog.GetFullPath()<<"\n";
	freopen(simLog.GetFullPath(), "w", stdout);

	//We can also drop in some penetrance files as well. But, for now, we'll be OK
}

AppController::SimParameters::SimParameters() : config("") {
	configuration = new AppConfig();
}

AppController::SimParameters::~SimParameters() {
	if (configuration)
		delete configuration;
}

AppController::~AppController()
{
	if (execLog)
		delete execLog;
	
	int gaCount = gaSettings.size();

	cout<<"Writing GA settings list to file\n";
	ofstream file("gaSettings.conf");
	for (int i=0; i<gaCount; i++) 
		file<<gaSettings[i];
	
}

void AppController::Refresh() {
	cout<<"Loading file: "<<config<<"\n";
	cout<<"Refresh()....need to load settings from the configuration\n";	
}


void AppController::SimParameters::DumpToFile() {
	cout<<"Saving configuration: "<<config<<"\n";
	cout<<"DumpToFile().....need to write settings to for configuration\n";
}

GrowthRate *AppController::SimParameters::GetGrowthRate() {
	return configuration->GetGrowthRate();
}

void AppController::SimParameters::SetGrowthRate(const char *growthCfg) {
	configuration->SetGrowthRate(growthCfg);
}


uint AppController::SimParameters::GetRandomSeed() {
	return configuration->GetRandomSeed();	
}

void AppController::SimParameters::SetRandomSeed(uint seed) {
	configuration->SetRandomSeed( seed );
}

uint AppController::SimParameters::GetDropCount() {
	return configuration->generalSettings.dropCount;
}

void AppController::SimParameters::SetDropCount(uint dc) {
	configuration->generalSettings.dropCount = dc;
}

uint AppController::SimParameters::GetDropFreq() {
	return configuration->generalSettings.dropFrequency;
}

void AppController::SimParameters::SetDropFreq(uint df) {
	configuration->generalSettings.dropFrequency = df;
}

uint AppController::SimParameters::GetDropInit() {
	return configuration->generalSettings.firstDropPoint;
}

void AppController::SimParameters::SetDropInit(uint di) {
	configuration->generalSettings.firstDropPoint = di;
}

void AppController::SimParameters::SetProjectName(const char *projectName) {
	configuration->generalSettings.outputName = projectName;
}

string AppController::SimParameters::GetProjectName() {
	return configuration->generalSettings.outputName;
}

uint AppController::SimParameters::GetSimultaneousChroms() {
	return PoolManager::simultaneousChrom;
}
	
void AppController::SimParameters::SetSimultaneousChroms(uint simChr){
	PoolManager::simultaneousChrom = simChr;
}
	
uint AppController::SimParameters::GetThreadsPerChrom() {
	return PoolManager::threadsPerChrom;
}

void AppController::SimParameters::SetThreadsPerChrom(uint threadCount) {
	PoolManager::threadsPerChrom = threadCount;
}

bool AppController::SimParameters::WriteConfiguration() {
	return configuration->WriteConfiguration(config.c_str());
}


bool AppController::SimParameters::WriteConfiguration(const char *config) {
	return configuration->WriteConfiguration(config);
}

void AppController::SimParameters::ClearLocusSelections() {
	configuration->generalSettings.locusSelections.clear();
}

void AppController::SimParameters::AddLocusSelector(LocusSelection &loc) {
	configuration->generalSettings.locusSelections[loc.GetLabel()] = loc;
}

LocusSelectionMap& AppController::SimParameters::GetLocusSelections() {
	PoolManager *pools = configuration->GetPoolManager();
	bool success = false;
	
	try {
		success = pools->InitializeLoci(configuration->generalSettings.firstGeneration, configuration->generalSettings.outputName.c_str(), configuration->generalSettings.doLoad);
	} catch (Utility::Exception::FileNotFound &e) { }

	if (success) {
		configuration->generalSettings.PartiallyReconcileSelections(*pools);
	}
	return configuration->generalSettings.locusSelections;
}



ChromPool *AppController::SimParameters::AddChromosomePool(uint blockCount, ChromPool::BlockDefinition& defBlock, const char *label) {
	return configuration->GetPoolManager()->AddChromosome(blockCount, defBlock, label);
}

ChromPool *AppController::SimParameters::LoadChromosomePool(const char *locSource, const char *label) {
	return configuration->GetPoolManager()->AddChromosome(locSource, label);
}

void AppController::SimParameters::DeleteChromosomePool(uint idx) {
	configuration->GetPoolManager()->DeleteChromosome(idx);
}

void AppController::SimParameters::SetDefaultBlock(ChromPool::BlockDefinition& defBlock) {
	configuration->SetDefaultBlock(defBlock);
}
ChromPool::BlockDefinition *AppController::SimParameters::GetDefaultBlock() {
	return configuration->GetDefaultBlock();
}

void AppController::SimParameters::SetFontFilename(const char *filename) {
	Simulation::Visualization::ImageParameters::font = filename;
}

string AppController::SimParameters::GetFontFilename() {
	return Simulation::Visualization::ImageParameters::font;
}

void AppController::SimParameters::SetMaxSnpDistance(uint distance) {
	Simulation::Visualization::LdPlotter::maxSnpDistance = distance;
}

uint AppController::SimParameters::GetMaxSnpDistance() {
	return Simulation::Visualization::LdPlotter::maxSnpDistance;
}

void AppController::SimParameters::SetStyleSheetFilename(const char *filename) {
	ChromPool::cssFilename = filename;
}

string AppController::SimParameters::GetStyleSheetFilename() {
	return ChromPool::cssFilename;
}

void AppController::SimParameters::SetReportBufferSize(uint buffSize) {
	ChromPool::highBufferSize = buffSize;
}

uint AppController::SimParameters::GetReportBufferSize() {
	return ChromPool::highBufferSize;
}

void AppController::SimParameters::SetTargetPopulation(uint targetPop) {
	ChromPool::targetPop = targetPop;
}

uint AppController::SimParameters::GetStartingGeneration() {
	return configuration->generalSettings.firstGeneration;
}

void AppController::SimParameters::SetStartingGeneration(uint gen) {
	configuration->generalSettings.firstGeneration = gen;
}

uint AppController::SimParameters::GetTargetPopulation() {
	return ChromPool::targetPop;
}

void AppController::SimParameters::WritePairwiseLD(bool doWrite) {
	ChromPool::writeLdTextReport = doWrite;
}

bool AppController::SimParameters::WritePairwiseLD() {
	return ChromPool::writeLdTextReport;
}

void AppController::SimParameters::PlotDPrime(bool doPlot) {
	ChromPool::writeDPrimePlots = doPlot;
}

bool AppController::SimParameters::PlotDPrime() {
	return ChromPool::writeDPrimePlots;
}

void AppController::SimParameters::PlotRSquared(bool doPlot) {
	ChromPool::writeRSquaredPlots = doPlot;
}	

bool AppController::SimParameters::PlotRSquared() {
	return ChromPool::writeRSquaredPlots;
}

void AppController::SimParameters::FastLD(bool isFast) {
	ChromPool::fastLD = isFast;
}

bool AppController::SimParameters::FastLD() {
	return ChromPool::fastLD;
}

void AppController::SummarizeSimulation(ostream &log) {
	//PoolManager *pools = parameters.configuration->GetPoolManager();
	parameters.configuration->GenerateReport(log);
}


PenetranceModel *AppController::SimParameters::GetDiseaseModel() {
	return configuration->statusSettings.model;
}

void AppController::SimParameters::SetDiseaseModel(PenetranceModel* model) {

	bool success = configuration->PostLoad();

	configuration->statusSettings.model = model;
}



void AppController::SimParameters::DefineDiseaseModel(const char *modelDetails) {
	configuration->statusSettings.DefineModel(modelDetails);
	
	try {
		configuration->ConfigureDiseaseModel(true);
	}
	catch (Utility::Exception::General& e) {
		wxMessageBox(_(e.GetErrorMessage().c_str()));
	}

}

void AppController::SummarizeDiseaseModel(ostream &os, vector<Locus*>&loci) {
	parameters.configuration->SummarizeDiseaseModel(os, loci);
}

void AppController::ClearLog() {
	execLog->ClearEntries();
}

void AppController::ClearCurrentLogEntries() {
	curEntries.clear();
}

size_t AppController::InitializeExecution(ostream &log) {
	simIsRunning = true;

	size_t startGen = parameters.GetStartingGeneration();
	if (!InitializeExecution(startGen, log)) {
		throw Utility::Exception::General("Unable to open project. Check that the files are available from this machine and that all disks are properly attached.");
	}
	simIsRunning = false;
	initializationComplete = true;
	return startGen;
}



bool AppController::LoadFileBasedChroms(vector<FileBasedChromosome *> &chroms) {
	assert(chroms.size() == 0);
	PoolManager *poolMgr = parameters.GetConfiguration()->GetPoolManager();
	PoolManager::Iterator itr = poolMgr->GetIterator();
	ChromPool *ch = itr.GetNext();
	const char *prj = parameters.GetProjectName().c_str();

	//Iterate over each of the chromosomes
	while (ch) {
		assert (ch->GetType() == ChromPool::FileBased);
		FileBasedChromosome *newChrom = new FileBasedChromosome(
				ch->GetID(), ch->GetLabel().c_str(), 
				ch->GenerateLocusFilename(prj).c_str(), 
				&(ch->GetLoci()));
		newChrom->Load();
		chroms.push_back(newChrom);
		ch = itr.GetNext();
	}
}


void AppController::InitExecLog() {
	if (execLog == NULL) {
		execLog = new ExecutionLog (parameters.configuration->GetAppConfigurationFilename().c_str());

		try {
			execLog->Initialize();
			if (pools)
				pools->CalculateAlleleFrequencies();
		} catch (Exception::FileIO& e) {
			wxMessageDialog dlg(NULL, _("Problems were encountered when trying to load information about previous runs associated with this configuration. As a result, not all data associated with this run will be available."));
			dlg.ShowModal();
		}
	}
}
bool AppController::InitializeExecution(uint startGeneration, ostream& log) {
	InitExecLog();

	if (!pools)
		pools = parameters.configuration->GetPoolManager();
	else
		pools->Close();

//	log<<"Initializing Pools\n";
	streambuf *oldCout = cout.rdbuf();
	cout.rdbuf(log.rdbuf());
	execLog->InitProject(parameters.configuration->generalSettings.seed, startGeneration, parameters.configuration->generalSettings.outputName.c_str());

//	if (startGeneration>0)
	parameters.configuration->generalSettings.firstGeneration = startGeneration;
	bool success = parameters.configuration->PostLoad();
	
	assert(pools->GetCurrentGeneration() == startGeneration);

	cout.rdbuf(oldCout);
//	log<<"\n\nChromosome Pools Initialized. \n";
	
	return success;
}



string AppController::OpenSummaryReport(ofstream &summary) {
	string reportFilename = "No_Report_Found";
	char *filename = new char[4096];
	simIsRunning = true;
	uint generation = parameters.configuration->GetPoolManager()->GetCurrentGeneration();
	string project = parameters.GetProjectName().c_str();
	string sampled = "";
	if (ChromPool::fastLD)
		sampled = "-sampled";
	sprintf(filename, "%s.Index%s-%d.ld.html", project.c_str(), sampled.c_str(), generation);
	reportFilename = filename;
	summary.open(filename);
	summary<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
	summary<<"<HTML><HEADER><TITLE> GenomeSIMLA Summary Report - "<<project<<" Generation "<<generation<<"</TITLE></HEADER>\n";
	summary<<"<BODY>\n";
	summary<<"<link rel=\"stylesheet\" type=\"text/css\" href='"<<ChromPool::cssFilename<<"'>\n<DIV class='image-frame'>\n";
	summary<<"<CENTER>\n<H1>GenomeSIMLA Summary Report</H1>\n";
	summary<<"\t<TABLE class='simple'><CAPTION>Execution Details</CAPTION>\n";
	summary<<"\t\t<TR><TH>Configuration</TH><TH>Project Name</TH><TH>Pool Size</TH></TR>\n";
	summary<<"\t\t<TR><TD>"<<parameters.GetConfigFilename()<<"</TD><TD>"<<project<<"</TD><TD>"<<parameters.GetGrowthRate()->GetPopulationSize(generation)<<"</TD></TR>\n\t</TABLE>\n";
	summary<<"\t<TABLE><CAPTION>Chromosomes</CAPTION>\n";

	delete[] filename;
	return reportFilename;
}

struct SimDetails {
	uint generationCount;
	AppController *controller;
	bool doLoadFirst;

	SimDetails(uint generationCount, AppController *ctrl, bool doLoadFirst) : 		
			generationCount(generationCount), controller(ctrl), doLoadFirst(doLoadFirst) {} 
	SimDetails() { exit(1);}
};
size_t AppController::RunSimulation(bool doLoadFirst) {
	size_t numGens = parameters.GetDropFreq();

	SimDetails *details = new SimDetails(numGens, this, doLoadFirst);
	pthread_create(&simThread, NULL, &ThreadedSimulation, (void*)details);

	return numGens;	
}
void AppController::PerformAnalysis() {
	pthread_create(&simThread, NULL, &ThreadedAnalysis, (void*)this);
}		
void AppController::WriteDatasets() {
	pthread_create(&simThread, NULL, &ThreadedDatasetGeneration, (void*)this);
}

void AppController::LockThread() {
	SIMLOCK;
}

void AppController::UnlockThead() {
	SIMUNLOCK;
}


void AppController::InitExecution() {
	pthread_create(&simThread, NULL, &ThreadedInitExecution, (void*)this);

}

void *AppController::ThreadedInitExecution(void *controller) {
	try {
		((AppController*)controller)->InitializeExecution(((AppController*)controller)->simOutput);
	} catch (Utility::Exception::General& e) {
		wxMessageBox(_(e.GetErrorMessage().c_str()));
		((AppController*)controller)->simIsRunning = false;
		((AppController*)controller)->initializationComplete = false;
		return NULL;
	}
	((AppController*)controller)->initializationComplete = true;
		
	return NULL;
}

void *AppController::ThreadedAnalysis(void *controller) {
	try {
	((AppController*)controller)->AnalyzeCurrentPool();
	} catch (Utility::Exception::General& e) {
		wxMessageBox(_(e.GetErrorMessage().c_str()));
		return NULL;
	}
	return NULL;
}


void *AppController::ThreadedDatasetGeneration(void *controller) {
	((AppController*)controller)->GenerateDatasets();
}



void AppController::GenerateDatasets() {
	simIsRunning = true;
	continueRunning = true;
	stringstream fileprefix;

	AppConfig* configuration = parameters.configuration;
	PoolManager *pools = configuration->GetPoolManager();


	size_t popSize = pools->GetExpressionCount();
	cout<<"Population Size: "<<popSize<<"\n";
	if (popSize < 35000) {
		wxMessageBox(wxT("The pool you have selected is too small to draw datasets from. Please select a pool with at least 35,000 chromosomes (we recommend 65,000)."), wxT("Pool Size"), wxICON_ERROR | wxOK);
		simIsRunning = false;
		return;
	}
	else if (popSize < 65000) 
		wxMessageBox(wxT("The pool you have selected is below the recommended size of 65,000. We highly recommend larger populations to ensure a reasonably diverse population."), wxT("Pool Size"), wxICON_WARNING | wxOK);

	PenetranceModel *model = configuration->statusSettings.model;
//	PenetranceModel *model = configuration->statusSettings.InitModel(pools, true);
//	model = configuration->statusSettings.LoadModel(configuration->generalSettings.outputName.c_str(), pools);
	

	fileprefix<<configuration->generalSettings.outputName<<"."<<pools->GetCurrentGeneration();
	
	uint sampleCount = configuration->datasetSettings.samples.size();
	uint datasetCount = configuration->datasetSettings.simSets;
	
	string sampleMarkerInfo;

	cout<<"Model use: "<<model->Details()<<"\n";
	percCompleted = 0.0;
	int totalDatasets = sampleCount * datasetCount;
	int datasetsCompleted = 0;
	//Iterate through the various samples that are to be created
	for (uint i=0; i<sampleCount && continueRunning; i++) {
		
		cout<<"\n\nWriting datasets for sample "<<configuration->datasetSettings.samples[i]->Details()<<"\n";
		for (uint p=0; p<datasetCount && continueRunning; p++) {
			cout<<"*";
			cout.flush();

			configuration->datasetSettings.samples[i]->BuildSample(*pools, model);
			configuration->datasetSettings.samples[i]->ApplyPhenocopyError(pools, model);			
			pools->ResolveGenotypes(configuration->generalSettings.outputName.c_str(), pools->GetCurrentGeneration(), configuration->datasetSettings.samples[i], model);
			configuration->datasetSettings.samples[i]->ApplyErrors(pools);
			configuration->datasetSettings.samples[i]->AppendGenotypeCountsToReport(model->GetDiseaseLoci());

			int snpCount = configuration->datasetSettings.samples[i]->GetSnpCount();
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
			filename<<configuration->datasetSettings.samples[i]->GetDescriptor();
			if (configuration->datasetSettings.binaryDatasets) {
				stringstream meta(filename.str().c_str());
				meta<<".meta";
				filename<<".bin";
				file.open(filename.str().c_str(), ios_base::out);
				ofstream metadata;
				if (configuration->datasetSettings.samples[i]->UsesMeta())
					metadata.open(meta.str().c_str(), ios_base::out);
				observations = configuration->datasetSettings.samples[i]->WriteBinaryDataset(metadata, file);
				file.close();
				metadata.close();
			} else {
				configuration->datasetSettings.samples[i]->WriteSnpDetails(filename.str().c_str(), pools);
				observations = configuration->datasetSettings.samples[i]->WriteDataset(filename.str().c_str(), gtCount);
			}
			percCompleted = ((float)++datasetsCompleted / (float)totalDatasets);
//			cout<<"--"<<percCompleted<<" ("<<datasetsCompleted<<"/"<<totalDatasets<<")\n";

			//configuration.datasetSettings.samples[i]->ReportGenotypeCounts(model->GetDiseaseLoci());
			configuration->datasetSettings.samples[i]->Purge();
		}

		cout<<"\n";

		SIMLOCK;
		configuration->datasetSettings.samples[i]->ReportGenotypeCounts(model->GetDiseaseLoci(), simOutput);
		SIMUNLOCK;
	}
	simIsRunning = false;

}

void AppController::CalculateAlleleFrequencies() {
	if (pools)
		pools->CalculateAlleleFrequencies();
}

void AppController::AnalyzeCurrentPool() {
	simIsRunning = true;
	string growthChart = parameters.configuration->GetGrowthChartFilename();
	if (PoolManager::generateOverviewLD) {

		if (curEntries.size() <= 0)
			return;
		string fastLD = "";

		if (parameters.FastLD())
			fastLD = "(Sampled)";
		else
			cout<<"Performing Complete LD Scan\n";
		SIMLOCK;
		simOutput<<"Writing Overview LD\n";
		SIMUNLOCK;

		ofstream summary;
		string mainReportFilename = OpenSummaryReport(summary);
		string locusReportFilename = pools->ProduceOverviewLD(parameters.GetProjectName().c_str(), false, true, summary, growthChart.c_str());
		ExecutionLog::LogEntry &current = curEntries[curEntries.size()-1];
		string lKey, iKey;

		current.reports[(string("Locus_Report")+fastLD)]=locusReportFilename;
		current.reports[(string("Index")+fastLD)]=mainReportFilename;
		current.sampledReport = parameters.FastLD();
		execLog->UpdateEntry(parameters.GetProjectName().c_str(), current);
		execLog->SaveFile();
	}
	simIsRunning = false;
}
void *AppController::ThreadedSimulation(void *simDetails) {
	SimDetails *details = (SimDetails*)simDetails;
	details->controller->RunSimulation(details->generationCount, details->doLoadFirst);
	return NULL;
}


/**
 * 
 * @param generationCount 
 * @param doLoadFirst 
 */
void AppController::RunSimulation(uint generationCount, bool doLoadFirst) {
	assert(pools);
	if (!initializationComplete) {
		return;
	}
	simIsRunning = true;
	uint curGen = pools->GetCurrentGeneration();
	uint curSize = pools->GetExpressionCount();
	uint firstDrop = parameters.GetDropInit();

	if (parameters.GetTargetPopulation() == 0 || parameters.GetTargetPopulation() > curSize) {
		cout<<"First Drop point: "<<firstDrop<<"\n";

		if (curGen < firstDrop) {
			generationCount = firstDrop - curGen;
			cout<<"Adjusting generation Count: "<<generationCount<<"\n";
		}
		SIMLOCK;

		simOutput<<"\nRunning Simulation from Generation "<<curGen<<" ("<<curSize<<") to "<<curGen+generationCount;
		SIMUNLOCK;
		pools->AdvanceGenerations(generationCount, parameters.configuration->GetGrowthRate(), doLoadFirst);
		int finalPop = pools->GetCurrentGeneration();

		if (finalPop != curGen+generationCount)
			simOutput<<"\nStopped at generation: "<<finalPop;
		execLog->StartRun(finalPop);

		SIMLOCK;
		string fastLD = "";
		if (parameters.FastLD())
			fastLD="(Sampled)";
		ExecutionLog::LogEntry entry = execLog->EndRun(finalPop, pools->GetExpressionCount(), "No_Report_Found", (string("Index")+fastLD).c_str());
		curEntries.push_back(entry);
		simOutput<<" ("<<pools->GetExpressionCount()<<")\n";
		
		SIMUNLOCK;
	}
	simIsRunning = false;
}


void AppController::AppendActiveEntry(ExecutionLog::LogEntry &entry) {
	curEntries.push_back(entry);
}
string AppController::GetSimOutput() {
	SIMLOCK;
	string buff = simOutput.str();
	simOutput.str("");

	assert(simOutput.str().length() == 0);
	SIMUNLOCK;
	return buff;
}
uint AppController::GetCurGeneration() {
	SIMLOCK;
	uint gen = pools->GetCurrentGeneration();
	SIMUNLOCK;
	return gen;
}

uint AppController::GetCurPopulation() {
	SIMLOCK;
	uint poolSize = pools->GetExpressionCount();
	SIMUNLOCK;
	return poolSize;
}


void AppController::HaltExecution() {
	ChromPool::continueRunning = continueRunning = false;
}

void AppController::GetSessionEntries(vector<ExecutionLog::LogEntry>& entries) {
	entries = curEntries;
}
bool AppController::GetProjectList(vector<string> &projects) {	
	InitExecLog();
	return execLog->GetProjectList(projects);
}

ExecutionLog::RunType *AppController::GetProjectEntries() {
	InitExecLog();
	return execLog->GetProjectEntries(parameters.GetProjectName().c_str());
}

ExecutionLog::RunType *AppController::GetProjectEntries(const char *project) {
	InitExecLog();
	return execLog->GetProjectEntries(project);
}

void AppController::UpdateCWD(const char *newCWD) {
	string oldCWD = wxGetCwd().c_str();

	PoolManager *poolMgr = parameters.GetConfiguration()->GetPoolManager();
	PoolManager::Iterator itr = poolMgr->GetIterator();
	ChromPool *ch = itr.GetNext();

	//Iterate over each of the chromosomes
	while (ch) {
		wxFileName relFilename (wxString(_(ch->GetLocSource().c_str())));	
		if (relFilename.MakeRelativeTo(newCWD))
			ch->SetLocSource(relFilename.GetFullPath().c_str());
		ch = itr.GetNext();
	}
}
void AppController::SetConfigurationFilename(const char *filename, bool doLoad) {
	if (execLog) {
		execLog->SaveFile();
		delete execLog;
	}
	execLog = NULL;
	wxFileName cfgFilename(filename);
	cfgFilename.MakeAbsolute();
	parameters.SetConfigFilename(cfgFilename.GetLongPath().c_str(), doLoad);
	UpdateCWD(cfgFilename.GetPath());
	wxFileName::SetCwd(cfgFilename.GetPath());


	delete execLog;
	execLog = NULL;
	InitExecLog();

}

}

}
