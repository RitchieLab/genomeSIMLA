//
// C++ Implementation: missingdataevaluation
//
// Description: 
//
//
// Author:  <>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "missingdataevaluation.h"
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

namespace MdrPDT {

namespace Evaluation {

MissingDataEvaluation::MissingDataEvaluation(float threshold) : threshold(threshold){	}


MissingDataEvaluation::~MissingDataEvaluation()
{
}

int MissingDataEvaluation::Analyze(std::ostream& os, GenotypeRepository& repo) {
	int foldCount = repo.GetFoldCount();
	int genotypeCount = repo.GetSnpCount();
	int indCount = repo.GetIndividualCount();
	int presentGenotypes = 0;

	int affected[genotypeCount];
	int unaffected[genotypeCount];
	vector<string> *labels = repo.GetSnpLabels();
	float dspCount = indCount/2.0;
	
	char *data=repo.GetSNPs();//+indCount;
	os<<"--------------- Missing Data -----------------\n"<<setw(10)<<left<<"SNP";
	for (int i=0; i<foldCount; i++) 
		os<<setw(12)<<right<<"Fold "<<setw(3)<<left<<i+1;
	os<<setw(15)<<right<<"Total Missing"<<"\n";
	for (int i=0;i<genotypeCount; i++) {
		char *folds=repo.GetFolds();
		int totalCount=0;
		memset(&affected, 0, foldCount*sizeof(int));
		memset(&unaffected, 0, foldCount*sizeof(int));
		for (int j=0; j<indCount; j+=2) {
			//Every other person is affected
			if ((*data++)>= 0) {
				affected[(int)*(folds++)]++;
				totalCount++;
			}
			if ((*data++)>=0){
				unaffected[(int)*(folds++)]++;
				totalCount++;
			}
		}
		float totalMissing = 0.0;
		if (totalCount > 0)
			totalMissing = (float)(2.0*dspCount - totalCount)/(2.0 * dspCount);
		if (totalMissing>threshold)
			repo.DoAnalyzeSNP(i+1, false);
		if (repo.DoAnalyzeSNP(i+1)) {
			os<<" ";
			presentGenotypes++;
		}
		else
			os<<"X";
		os<<setw(9)<<left<<(*labels)[i];
		for (int f=0; f<foldCount; f++) {
			stringstream ss;
			ss<<affected[f]<<"A/"<<unaffected[f]<<"U";
			//ss<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(2)<<( ( (dspCount - (float)affected[f]) / dspCount ) * 100.0)
			//	<<"/"<<setiosflags(ios::fixed|ios::showpoint)<<(((dspCount - (float)unaffected[f]) / dspCount)*100.0)<<"% ";
			os<<right<<setw(15)<<ss.str();
		}
		os<<"      "<<setprecision(4)<<(totalMissing*100.0)<<"%\n";
	}
	return presentGenotypes;
}

}

}
