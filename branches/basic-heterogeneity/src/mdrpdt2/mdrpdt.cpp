//
// C++ Implementation: mdrpdt
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <string>
#include <stdlib.h>
#include "tstatistic.h"

#include "familyrepository.h"
#include "ptestdistribution.h"
#include "appmdrpdt.h"

using namespace std;
using namespace MdrPDT;

#include "utility/random.h"

int main(int argc, char *argv[])	{
	string cfgFilename;
	int reportSize = 5;

	Evaluation::MatchedOddsRatio::Verbose = false;
	string dataset;

	AppMdrPDT app;					///<The application object
	if (!app.ParseCmdLine(argc, argv))
		exit(1);
	
	ApplicationConfiguration *cfg = ConfigurationReader::appConfiguration;
	reportSize = cfg->GetInteger("REPORTMODELCOUNT");

	//The report which keeps all of our results from the real run
	Evaluation::TFinalReport topModels(Evaluation::EvaluationMethod::maxModelSize, app.GetFoldCount(), reportSize);
	app.LoadData(cfg->GetLine("INPUTFILE").c_str());
	app.EvaluateDataset(topModels);
	
	if (app.IsMasterNode() ){	
		//Let's print out the details for the top model(s) 
		topModels.GenerateDetailedReport(0, cout, app.GetEvaluationMethod());
	}

	int ptestCount = cfg->GetInteger("PTEST_COUNT");
	if (app.IsMasterNode()) {
		cout<<"\n\n-------------- Running PTests ("<<ptestCount<<") ------------\n";
		//Set up the best score so that we can short circuit the ptests
		string bestModel;
		float bestMOR = topModels.GetBestMOR(bestModel);
		app.SetBestModel(bestMOR, bestModel.c_str(), ptestCount);
	}

	//Now, we'll run our ptests
	Distribution::PTestDistribution *dist = app.RunPTests(ptestCount);

	if (app.IsMasterNode()) {
		string distFilename = app.GetFilename(cfg->GetLine("EXT_DISTRIBUTION").c_str());
		ofstream distFile(distFilename.c_str());
		dist->ReportDistribution(0,distFile);
		distFile.close();
		cout<<"\n\nDistribution Size:	"<<dist->GetDistributionSize()<<"\n";
		//Produce the summary report
		topModels.GenerateSummaryReport(cout, dist);
	}
	if (dist)
		delete dist;
  	return EXIT_SUCCESS;
}
