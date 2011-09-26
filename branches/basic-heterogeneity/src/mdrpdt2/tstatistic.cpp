//
// C++ Implementation: tstastic
//
// Description: 
//
//
// Author:  <Eric Torstenson>, (C) Marylyn Ritchie 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tstatistic.h"
#include <sstream>
#include "matchedoddsratio.h"
#include "foldstatistic.h"

using namespace std;

namespace MdrPDT {

using namespace Distribution;


namespace Evaluation { 

int EvaluationMethod::maxModelSize = 3;
int EvaluationMethod::minModelSize = 1;

TStatistic::TStatistic(int xvCount, int pedCount, int indCount) : 
		xvCount(xvCount), pedCount(pedCount), indCount(indCount), snpCount(0), folds(NULL), pedIDs(NULL) {  }

TStatistic::~TStatistic() { }


void TStatistic::BasicEval(GenotypeRepository *rep, TFinalReport& report, vector<char *> snps, vector<int> modelIDs, int snpID) {
	int modelSize = snps.size();
	if (modelSize > maxModelSize || snpID > snpCount) 
		return;
		
	modelSize++;
	while (snpID <= snpCount) {
		if (rep->DoAnalyzeSNP(snpID)) {
			snps.push_back(rep->GetSNP(snpID-1));
			modelIDs.push_back(snpID);
	
			if (modelSize >= minModelSize) {
				PdtModel curModel;
				curModel.SetSnpLabels(rep->GetSnpLabels());
				curModel.SetDSPCounts(rep->GetDSPCounts());

				for (int i=0; i<modelSize; i++) 
					curModel.AddSnp(modelIDs[i], snps[i]);
				curModel.InitPedigree(indCount, pedCount, pedIDs);
	
	
				curModel.CrossValidations(xvCount, folds);			
				EvaluateModel(curModel);
				report.AddModelResult(curModel);
			}
		
			//Recursively finish off this line
			if (modelSize < maxModelSize) 
				BasicEval(rep, report, snps, modelIDs, snpID+1);

			modelIDs.pop_back();
			snps.pop_back();
		}
		snpID++;
	}	
}

void TStatistic::BasicEval(GenotypeRepository *rep, TFinalReport& report) {
	snpCount = rep->GetSnpCount();
	folds = rep->GetFolds();
	pedIDs = rep->GetPedigreeIDs();

	vector<char *> snps;
	vector<int> modelIDs;
	BasicEval(rep, report, snps, modelIDs, 1);
}

void TStatistic::EvaluateModel(GenotypeRepository *rep, const char *model) {
	snpCount = rep->GetSnpCount();
	folds = rep->GetFolds();
	pedIDs = rep->GetPedigreeIDs();

	string sModel = model;
	boost::char_separator<char> sep("x");
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	tokenizer snps(sModel, sep);
	tokenizer::iterator itr = snps.begin();
	tokenizer::iterator end = snps.end();
	
	PdtModel curModel;
	curModel.SetSnpLabels(rep->GetSnpLabels());
	curModel.SetDSPCounts(rep->GetDSPCounts());
	while (itr!=end) {
		int snpID = atoi(itr->c_str());
		if (snpID < 1) {
			cerr<<"Unable to work with model, "<<model<<"\n";
			exit(1);
		}
		curModel.AddSnp(snpID, rep->GetSNP(snpID-1));
		itr++;
	}

	curModel.InitPedigree(indCount, pedCount, pedIDs);
	curModel.CrossValidations(xvCount, folds);		

	//We have to evaluate the model before we can describe it...counts are not done in the verbose version
	EvaluateModel(curModel);

	if (xvCount > 1) {
		for (int i=0; i<xvCount; i++) {
			EvaluateModelVerbose(curModel, i, cout);
		}
	}
	else {
		EvaluateModelVerbose(curModel, cout);
	}
	exit(0);
}

/**
 * @brief returns genotype at locus, locus for the multilocus genotype, genotype
 */
void TStatistic::GetGenotype(int genotype, int locusCount, deque<int>& genotypes) {
	genotypes.clear();
	int multiplier = (int)pow(3, locusCount - 1);
	//int multiplier = 1;
	int remainder = genotype;
	for (int i=0;i<locusCount; i++) {
		int genotype = remainder/multiplier;
		remainder-=(genotype*multiplier);
		genotypes.push_front(genotype);
		multiplier/=3;
	}
	
}



void TStatistic::EvaluateModel(PdtModel& report) {
	report.InitModelCounts();
	int foldCount			= report.FoldCount();			
	
	OverallCounter& genotypeCounter = report.GetGenotypeCounter();
	
	
	if (foldCount > 1) {
		for (int i=0; i<foldCount; i++) {
			tcalc training;
			tcalc testing;
			genotypeCounter.EvaluateFold(i, training, testing);
			//genotypeCounter.ReportFold(i, training, testing);
			report.AddFold(training,testing);			
		}
	}
	else {
		vector<int> hrCells;
		tcalc results = genotypeCounter.EvaluateAsOne(hrCells);
		report.AddFold(results, results);
	}
}

void TStatistic::EvaluateModelVerbose(PdtModel& report, int fold, ostream& os) {
	OverallCounter& genotypeCounter = report.GetGenotypeCounter();
	
	//genotypeCounter.DebugReport(cout);

	int genotypeCount = genotypeCounter.genotypeCount;
	vector<int> hrCells;
	///<The snp IDs of interest
	vector<int> snpIDs = report.GetModelIDs();			
	//First, we need to gather the details we are about to use
	cout<<"\n\nFold "<<fold+1<<" Details ( "<<report.GetModelID()<<" ):\n";
	int totalAffected	= 0;
	int totalUnaffected	= 0;
	genotypeCounter.GetTotalGenotypeCounts(totalAffected, totalUnaffected);

	stringstream gtLabel;
	vector<int> labelWidths;
	int totalLabelWidth = 0;
	for (int n=0; n<report.ModelSize(); n++) {
		string label = report.GetSnpLabel(n);
		int w = label.length()+2;
		totalLabelWidth+=w;
		labelWidths.push_back(w);
		gtLabel<<setw(w)<<right<<report.GetSnpLabel(n)<<" ";
	}
	if (totalLabelWidth < 9)
		totalLabelWidth=9;

	os<<"\n"<<setw(totalLabelWidth)<<" Genotype ";
	os<<"  "<<setw(12)<<"Affected"<<setw(12)<<"Unaffected"<<"\n";

	os<<setw(totalLabelWidth)<<right<<gtLabel.str()<<" "<<setw(12)<<right<<"Trn/Test"<<" "<<setw(12)<<right<<"Trn/Test";
	os<<"  "<<setw(12)<<"Total";
	os<<"  "<<setw(12)<<"Ratio"<<"\n";

	deque<int> snpGenotypes;

	string tag;
	tcalc training;
	tcalc testing;	

	int affectedCount, unaffectedCount;
	genotypeCounter.GetTotalGenotypeCounts(fold, affectedCount, unaffectedCount);

	cout<<"  ------------------------------------------------------------------------\n";
	hrCells = genotypeCounter.EvaluateFold(fold, training, testing);
	for (int g=0; g<genotypeCount; g++) {
		GetGenotype(g, report.ModelSize(), snpGenotypes);

		if (find(hrCells.begin(), hrCells.end(), g) != hrCells.end())
			tag="*";
		else
			tag=" ";
		
		int trainingAff, trainingUnaff, testingAff, testingUnaff;
		genotypeCounter.GetGenotypeCounts(fold, g, trainingAff, 
				trainingUnaff, testingAff, testingUnaff);

		stringstream genoID;
		stringstream affValue;
		stringstream unaffValue;
		for (int i=0;i<report.ModelSize(); i++) 
			genoID<<setw(labelWidths[i])<<right<<report.GetGenotypeEncoding(i, snpGenotypes[i]);
		genoID<<tag;
		affValue<<trainingAff<<"/"<<testingAff;
		unaffValue<<trainingUnaff<<"/"<<testingUnaff;
		os<<setw(totalLabelWidth)<<genoID.str()<<" "
			<<setw(12)<<affValue.str()<<" "
			<<setw(12)<<unaffValue.str()<<" "
			<<setw(12)<<(trainingAff+trainingUnaff);
		os<<setw(12)<<setiosflags(ios::fixed | ios::showpoint)
			 <<setprecision(4)<<((float)trainingAff/(float)trainingUnaff)<<endl;
	}
	float dspCount = report.IndividualCount(fold)/200.0;
	os<<"  ------------------------------------------------------------------------\n";
	os<<setw(totalLabelWidth)<<" "<<"  "<<setw(12)<<affectedCount<<" "<<setw(12)<<unaffectedCount<<" "
		<<setw(12)<<affectedCount+unaffectedCount<<"\n";
	os<<setw(totalLabelWidth)<<"Missing: "<<"  "
		<<right<<setw(12)<<setprecision(2)<<100.0-((float)affectedCount/dspCount)<<"%"
		<<right<<setw(12)<<setprecision(2)<<100.0-((float)unaffectedCount/dspCount)<<"%\n\n";

	os<<setw(35)<<right<<"T-Statistic (Training) "<<training.GetTStatistic()<<"\n";
	os<<setw(35)<<right<<"T-Statistic (Testing) "<<testing.GetTStatistic()<<"\n";
	MatchedOddsRatio mor = report.EvaluateMOR(fold);
	os<<setw(35)<<right<<"Matched Odds Ratio "<<mor.GetRatio()<<"\n";
	//This is a bit wrong, but, it's really the only way to get the MOR inside the report
	
	
}

void TStatistic::EvaluateModelVerbose(PdtModel& report, ostream& os) {
	OverallCounter& genotypeCounter = report.GetGenotypeCounter();
	
	//genotypeCounter.DebugReport(cout);

	int genotypeCount = genotypeCounter.genotypeCount;
	vector<int> hrCells;
	///<The snp IDs of interest
	vector<int> snpIDs = report.GetModelIDs();			
	//First, we need to gather the details we are about to use
	cout<<"\n\n Details ( "<<report.GetModelID()<<" ):\n";
	int totalAffected	= 0;
	int totalUnaffected	= 0;
	genotypeCounter.GetTotalGenotypeCounts(totalAffected, totalUnaffected);

	stringstream gtLabel;
	vector<int> labelWidths;
	int totalLabelWidth = 0;
	for (int n=0; n<report.ModelSize(); n++) {
		string label = report.GetSnpLabel(n);
		int w = label.length()+2;
		totalLabelWidth+=w;
		labelWidths.push_back(w);
		gtLabel<<setw(w)<<right<<report.GetSnpLabel(n)<<" ";
	}
	if (totalLabelWidth < 9)
		totalLabelWidth=9;
	os<<"\n"<<setw(totalLabelWidth)<<right<<" Genotype ";
	os<<"  "<<setw(12)<<" "<<setw(12)<<" "<<"\n";

	os<<setw(totalLabelWidth)<<right<<gtLabel.str()<<" "<<setw(12)<<right<<"Affected"<<" "<<setw(12)<<right<<"Unaffected";
	os<<"  "<<setw(12)<<"Total";
	os<<"  "<<setw(12)<<"Ratio"<<"\n";

	deque<int> snpGenotypes;

	string tag;
	tcalc tStat;

	int affectedCount, unaffectedCount;
	genotypeCounter.GetTotalGenotypeCounts(affectedCount, unaffectedCount);

	cout<<"  ------------------------------------------------------------------------\n";
	tStat = genotypeCounter.EvaluateAsOne(hrCells);
	for (int g=0; g<genotypeCount; g++) {
		GetGenotype(g, report.ModelSize(), snpGenotypes);

		if (find(hrCells.begin(), hrCells.end(), g) != hrCells.end())
			tag="*";
		else
			tag=" ";
		
		int aff, unaff;
		genotypeCounter.GetGenotypeCounts(g, aff, unaff);

		stringstream genoID;
		for (int i=0;i<report.ModelSize(); i++) 
			genoID<<setw(labelWidths[i])<<right<<report.GetGenotypeEncoding(i, snpGenotypes[i]);
		genoID<<tag;

		os<<setw(totalLabelWidth)<<genoID.str()<<" "
			<<setw(12)<<aff<<" "
			<<setw(12)<<unaff<<" "
			<<setw(12)<<(aff+unaff);
		os<<setw(12)<<setiosflags(ios::fixed | ios::showpoint)
			 <<setprecision(4)<<((float)aff/(float)unaff)<<endl;
	}
	float dspCount = report.IndividualCount()/200.0;
	os<<"  ------------------------------------------------------------------------\n";
	os<<setw(totalLabelWidth)<<" "<<"  "<<setw(12)<<affectedCount<<" "<<setw(12)<<unaffectedCount<<" "
		<<setw(12)<<affectedCount+unaffectedCount<<"\n";
	os<<setw(totalLabelWidth)<<"Missing: "<<"  "
		<<right<<setw(12)<<setprecision(2)<<((float)affectedCount/dspCount)<<"%"
		<<right<<setw(12)<<setprecision(2)<<((float)unaffectedCount/dspCount)<<"%\n\n";
	os<<setw(35)<<right<<"T-Statistic "<<tStat.GetTStatistic()<<"\n";
	MatchedOddsRatio mor = report.EvaluateMOR();
	os<<setw(35)<<right<<"Matched Odds Ratio "<<mor.GetRatio()<<"\n";
	//This is a bit wrong, but, it's really the only way to get the MOR inside the report
	
	
}



}

}



