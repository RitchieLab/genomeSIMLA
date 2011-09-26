/*
 *  genotypecounter.h
 *  Project
 *
 *  Created by Eric Torstenson on 8/27/08.
 *  Copyright 2008 Marylyn Ritchie. 
 *
 *  See COPYING file that comes with this distribution
 *
 */
 #ifndef GENOTYPECOUNTER_H
 #define GENOTYPECOUNTER_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>


#include "matchedoddsratio.h"

namespace MdrPDT {
namespace Evaluation {

using namespace std;

/**
 * @brief This represents the score for a given fold (or the entire set) for a specific model
 */
struct tcalc {
	tcalc() : sumD(0), sumD2(0) { }

	void AppendD(int d) {
		sumD+=d;
		sumD2+=(d*d);
	}
	
	tcalc& operator+=(const tcalc& other) {
		sumD+=other.sumD;
		sumD2+=other.sumD2;
		return *this;
	}
	
	void Reset() { sumD=0; sumD2=0; }
	
	float GetTStatistic() const {
		return (float)sumD/(sqrt(sumD2));
	}
	
	void Report(ostream& os) {
		os<<setw(8)<<sumD<<setw(8)<<sumD2<<setw(8)<<GetTStatistic();
	}
	int sumD;
	int sumD2;
};


struct GenotypeConverter {
	GenotypeConverter() : indCount(0) { }
	GenotypeConverter(int modelSize, int indCount) : indCount(indCount) {
		genotypeCount=3;
		multiplier.push_back(1);
		for (int i=1; i<modelSize; i++) {
			multiplier.push_back(genotypeCount);
			genotypeCount*=3;
		}
	}
	bool GetGenotype(vector<char *>genotypes, int offset, int &aff, int &unaff) {
		aff=unaff=0;
		int modelSize = multiplier.size();
		for (int i=0; i<modelSize; i++) {
			int gtA = genotypes[i][offset]-1;
			int gtU = genotypes[i][offset + 1] - 1;

			if (gtA >= 0 && gtU >=0) {
				aff += gtA * multiplier[i];
				unaff += gtU * multiplier[i];
			}
			else {
				aff = unaff = -1;
				break;
			}
		}
		return aff >= 0 && unaff >=0;
	}

	int GetModelSize() { return multiplier.size(); }
	vector<int> multiplier;		
	int genotypeCount;
	int indCount;
};
/** 
 * @brief Keeps up with the number of observations for each genotype
 */
class GenotypeCounts {
public:
	GenotypeCounts(int genotypeCount) : sum(0), genotypeCounts(genotypeCount, 0) { }
	int GetCount(int idx) { return genotypeCounts[idx]; }
	int GetCellCount() { return genotypeCounts.size(); }
	int GetTotal() { return sum; }
	void AddCount(int genotype) { genotypeCounts[genotype]++;}
	void Report(ostream& os) {
		for (size_t i=0; i<genotypeCounts.size(); i++) {
			os<<" "<<genotypeCounts[i]<<" ";
		}
	}
protected:
	int sum;								///<Count of all genotypes observed(individuals)
	vector<int> genotypeCounts;				///<Sum for each genotype
};

/**
 * @brief Abstract out the concept of counting the different genotypes for all the different striations
 */
class PedigreeCounter {
public:
	PedigreeCounter(int genotypeCount) : sum(0), affectedGenotypes(genotypeCount), unaffectedGenotypes(genotypeCount) { }
	/**
	 * @brief Increments the counts for each of the specified genotypes
	 */
	void AddCounts(int aff, int unaff) { 
		affectedGenotypes.AddCount(aff);
		unaffectedGenotypes.AddCount(unaff);
		sum+=2;
	}
		
	/**
	 * @brief Returns the number of individuals with the given genotype (affCount, unaffCount and returns difference)
	 */
	int EvaluateGenotype(int genotype, int &affCount, int &unaffCount) {
		affCount=affectedGenotypes.GetCount(genotype);
		unaffCount = unaffectedGenotypes.GetCount(genotype);
		return affCount-unaffCount;
	}
	void EvaluateMOR(MatchedOddsRatio &mor, vector<int>& hrCells) {
		int truePos=0, falsePos=0, trueNeg=0, falseNeg=0;
		int genotypeCount = affectedGenotypes.GetCellCount();
		for (int i=0; i<genotypeCount; i++) {
			//HR Cells
			if (find(hrCells.begin(), hrCells.end(), i) != hrCells.end()) {
				truePos+=affectedGenotypes.GetCount(i);
				falsePos+=unaffectedGenotypes.GetCount(i);
			}else {
				falseNeg+=affectedGenotypes.GetCount(i);
				trueNeg+=unaffectedGenotypes.GetCount(i);
			}
		}
		mor.AddFamily(truePos, falsePos, falseNeg, trueNeg);
		
	}
	
	void GetGenotypeCounts(int genotype, int& aff, int& unaff) {
		aff=affectedGenotypes.GetCount(genotype);
		unaff=unaffectedGenotypes.GetCount(genotype);
	}

	void Report(ostream& os) {
		affectedGenotypes.Report(os);
		os<<" | ";
		unaffectedGenotypes.Report(os);
	}
	
	int Sum() { return sum; }
protected:
	int sum;								///<Observed for whole pedigree
	GenotypeCounts affectedGenotypes;		///<Counts for affected individuals at each genotype
	GenotypeCounts unaffectedGenotypes;		///<Counts for unaffected individuals at each genotype
};

class GenotypeCounter {
public:
	GenotypeCounter(int pedigreeCount, int genotypeCount) 
		: pedigreeCounters(pedigreeCount, genotypeCount), totalAffecteds(0), totalUnaffecteds(0) { }

	void AddCounts(int pedID, int affGenotype, int unaffGenotype) {
		pedigreeCounters[pedID].AddCounts(affGenotype, unaffGenotype);
		totalAffecteds++;
		totalUnaffecteds++;
	}
	void GetGenotypeCounts(int& aff, int& unaff) {
		aff = totalAffecteds;
		unaff = totalUnaffecteds;
	}
	void GetGenotypeCounts(int genotype, int& aff, int& unaff) {
		vector<PedigreeCounter>::iterator itr = pedigreeCounters.begin();
		vector<PedigreeCounter>::iterator end = pedigreeCounters.end();
		aff=unaff=0;
		while (itr!=end) {
			int a, u;
			itr->GetGenotypeCounts(genotype, a, u);
			aff+=a;
			unaff+=u;
			itr++;
		}
	}

	//Iterate overall Pedigrees and sum up their D statistics
	tcalc EvaluateGenotype(vector<int> hrCells) {
		int aff, unaff;
		tcalc sumD;
		vector<PedigreeCounter>::iterator itr = pedigreeCounters.begin();
		vector<PedigreeCounter>::iterator end = pedigreeCounters.end();

		while (itr!=end) {
			//Work through each part of the model
			tcalc localSumD;
			vector<int>::iterator gt=hrCells.begin();
			vector<int>::iterator endGT=hrCells.end();

			aff=0;
			unaff=0;
			while (gt!=endGT) {
				int lA, lU;
				itr->EvaluateGenotype(*gt, lA, lU);
				aff+=lA;
				unaff+=lU;
				gt++;
			}
			itr++;
			sumD.AppendD(aff-unaff);
		}
		return sumD;
	}
	
	MatchedOddsRatio EvaluateMOR(vector<int>& hrCells) {
		MatchedOddsRatio mor;
		vector<PedigreeCounter>::iterator itr = pedigreeCounters.begin();
		vector<PedigreeCounter>::iterator end = pedigreeCounters.end();
		while (itr!=end) {
			itr->EvaluateMOR(mor, hrCells);
			itr++;
		}
		return mor;
	}
	
	void Report(ostream& os) {
		for (size_t i=0; i<pedigreeCounters.size(); i++) {
			os<<i<<"\t";
			pedigreeCounters[i].Report(os);
			os<<"\n";
		}
	}
protected:
	vector<PedigreeCounter> pedigreeCounters;
	int totalAffecteds;
	int totalUnaffecteds;
};

class OverallCounter {
public:
	OverallCounter() : genotypeCount(0), overallAffected(0), overallUnaffected(0) { }
	OverallCounter(int foldCount, int pedigreeCount, int genotypeCount) : genotypeCount(genotypeCount),
			totalAffected(genotypeCount, 0), totalUnaffected(genotypeCount, 0), overallAffected(0), overallUnaffected(0), folds(foldCount, GenotypeCounter(pedigreeCount, genotypeCount)) {}

	void AddCounts(int foldID, int pedID, int gtAff, int gtUnaff) {
		totalAffected[gtAff]++;
		totalUnaffected[gtUnaff]++;
		overallAffected++;
		overallUnaffected++;
		folds[foldID].AddCounts(pedID, gtAff, gtUnaff);
	}
	tcalc EvaluateAsOne(vector<int>& hrCells) {
		HrCells(hrCells);
		tcalc result;							///<For non-CV, there is no training/testing
		
		vector<GenotypeCounter>::iterator itr = folds.begin();
		vector<GenotypeCounter>::iterator end = folds.end();
		//Sum up all results from each fold
		while (itr != end) {
			result+=itr->EvaluateGenotype(hrCells);
			itr++;
		}
//cout<<"--->"<<hrCells.size()<<" "<<result.sumD<<"  "<<result.sumD2<<"\n";
		return result;
	}

	MatchedOddsRatio EvaluateMOR() {
		vector<int>hrCells;		
		tcalc tScore = EvaluateAsOne(hrCells);
		MatchedOddsRatio mor;

		for (uint i=0; i<folds.size(); i++) {
			mor = mor + folds[i].EvaluateMOR(hrCells);
		}
		return mor;
	}

	MatchedOddsRatio EvaluateMOR(int fold) {
		if (folds.size() == 1)
			return EvaluateMOR();

		tcalc trainingT;
		tcalc testingT;
		vector<int>hrCells = EvaluateFold(fold, trainingT, testingT);
		MatchedOddsRatio mor = folds[fold].EvaluateMOR(hrCells);
		return mor;
	}
	
	MatchedOddsRatio EvaluateBiasedMOR(int fold) {
		tcalc trainingT;
		tcalc testingT;
		vector<int>hrCells = EvaluateFold(fold, trainingT, testingT);
		MatchedOddsRatio mor=folds[fold].EvaluateMOR(hrCells);
		return mor;
	}

	void HrCells(int fold, vector<int>& hrCells) {
		int aff, unaff;
		hrCells.clear();
		for (int gt=0; gt<genotypeCount; gt++) {
			
			folds[fold].GetGenotypeCounts(gt, aff, unaff);
			if (totalAffected[gt]-aff >= totalUnaffected[gt]-unaff)
				hrCells.push_back(gt);
		}
	}

	void HrCells(vector<int>& hrCells) {
		int foldCount = folds.size();
		hrCells.clear();
		for (int gt=0; gt<genotypeCount; gt++) {
			for (int f=0; f<foldCount; f++) {
				int lAff=0, lUnaff=0;

				folds[f].GetGenotypeCounts(gt, lAff, lUnaff);
				if (lAff>=lUnaff)
					hrCells.push_back(gt);
			}
		}	
		
	}
	/**
	 * Returns the vector of HR cells
	 */
	vector<int> EvaluateFold(int fold, tcalc& training, tcalc& testing) {
		training.Reset();
		testing.Reset();
		vector<int> hrCells;
		HrCells(fold, hrCells);
		int foldCount = folds.size();
			
		for (int f=0; f<foldCount; f++) {
			if (f==fold)
				testing=folds[f].EvaluateGenotype(hrCells);
			else	
				training+=folds[f].EvaluateGenotype(hrCells);
		}
		return hrCells;
	}
	
	void DebugReport(ostream& os) {
		tcalc trainingT;
		tcalc testingT;
		vector<int> hrCells;
		for (size_t i=0; i<folds.size();i++) {
			os<<"\n\nFold ("<<i<<"):\n\tAffected  | Unaffected\n\t";

			hrCells = EvaluateFold(i, trainingT, testingT);
			for (int a=0; a<2; a++) {
				for (int g=0; g<genotypeCount; g++) {
					if (find(hrCells.begin(), hrCells.end(), g) != hrCells.end())
						os<<"*";
					else
						os<<" ";
					os<<i<<" ";
				}
				os<<" | ";
			}
			os<<"\n";
			folds[i].Report(os);
			os<<"----------------------------------\n";
			os<<"Training: "<<trainingT.sumD<<"\t"<<trainingT.sumD2<<"\t"<<trainingT.GetTStatistic()<<"\n";
			os<<"Testing:  "<<testingT.sumD<<"\t"<<testingT.sumD2<<"\t"<<testingT.GetTStatistic()<<"\n";
		}
	}

	void GetTotalGenotypeCounts(int fold, int& aff, int &unaff) {
		return folds[fold].GetGenotypeCounts(aff, unaff);
	}

	void GetTotalGenotypeCounts(int& aff, int& unaff) {
		aff=overallAffected;
		unaff=overallUnaffected;
	}

	void GetGenotypeCounts(int genotype, int& aff, int& unaff) {
		aff=totalAffected[genotype];
		unaff=totalUnaffected[genotype];
	}

	void GetGenotypeCounts(int fold, int genotype, int& trAff, 
			int& trUnaff, int& tsAff, int& tsUnaff) {
		trAff=trUnaff=tsAff=tsUnaff=0;
		int foldCount = folds.size();
		for (int f=0; f<foldCount; f++) {
			if (f==fold) 
				folds[f].GetGenotypeCounts(genotype, tsAff, tsUnaff);
			else {
				int a, u;
				folds[f].GetGenotypeCounts(genotype, a, u);
				trAff+=a;
				trUnaff+=u;
			}
		}
	}

	/**
	 * @brief Returns the counts of the testing set for a given genotype at the given fold
	 */
	void GetGenotypeCounts(int fold, int genotype, int& aff, int& unaff) {
		folds[fold].GetGenotypeCounts(genotype, aff, unaff);
	}
	
	void ReportFold(int fold, tcalc& training, tcalc& testing) {
		ostream &os = cout;
		os<<fold<<"\t"<<training.GetTStatistic()<<"\t"<<testing.GetTStatistic();
	}
	int genotypeCount;
	vector<int> totalAffected;
	vector<int> totalUnaffected;
	int overallAffected;
	int overallUnaffected;
	vector<GenotypeCounter> folds;
};
}
}
#endif
