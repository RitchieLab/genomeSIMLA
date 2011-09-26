#ifndef DISEASELOCI
#define DISEASELOCI

#include <math.h>
#include <vector>
#include <iostream>

namespace Simulation {
namespace StatusModel {

using namespace std;


struct DiseaseLocus {
	int chromosome;					///<The chromosome the locus is located on
	int locusIdx;					///<The index of the locus
	float alFreq1;					///<Allele frequency
	float alFreq2;					///<Allele frequency
	string label;					///<Textual identification
	
	DiseaseLocus() : chromosome(-2), locusIdx(-2), alFreq1(0.0), alFreq2(0.0) { }
	DiseaseLocus(const char *label, int c, int i, float f1, float f2) : chromosome(c), locusIdx(i), alFreq1(f1), alFreq2(f2), label(label) { }
	DiseaseLocus(const DiseaseLocus& other) : 
							chromosome(other.chromosome),
							locusIdx(other.locusIdx),
							alFreq1(other.alFreq1),
							alFreq2(other.alFreq2),
							label(other.label) { 
		//cout<<"DiseaseLocus(dl):\n";
		//cout<<"\tChrom: "<<other.chromosome<<"\n\tLoc:  "<<other.locusIdx<<"\n\t"<<other.alFreq1<<"\n";
	}

	bool operator<(const DiseaseLocus& other) const {
		if (chromosome == other.chromosome) {
			return locusIdx < other.locusIdx;
		}
		return chromosome < other.chromosome;
	}

	bool operator==(const DiseaseLocus& other) const {
		if (chromosome == other.chromosome)
			return locusIdx == other.locusIdx;
		return false;
	}
	/**
	 * @Brief returns true of contents <I>appear</I> to be legitimate. 
	 */
#ifdef USE_XY
	bool ValidIdx() { return (chromosome > -2) && locusIdx > -1; }
#else
	bool ValidIdx() { return chromosome > -1 && locusIdx > -1; }
#endif //USE_XY
	/**
	 * @brief Checks that the allele frequencies fit together well enough
	 */
	bool IsValid()  { return ValidIdx() && alFreq1 + alFreq2 > 0.9999 && alFreq1 + alFreq2 < 1.0001; }

	/**
	 * @brief Returns that the frequencies are valid (sum to 1.0)
	 */
	bool Verify() {  return fabs(1.0 - (alFreq1 + alFreq2)) < 0.001; }

	void Report() { cout<<"\t"<<chromosome<<":"<<locusIdx<<" "<<alFreq1<<" "<<alFreq2<<" "<<fabs(1.0 - (alFreq1 + alFreq2))<<"\n"; }

	/**
	 * @Brief Returns the minor allele frequency
	 */
	float GetMinAlleleFreq() { if (alFreq1 < alFreq2) return alFreq1; else return alFreq2; }

};


typedef std::vector<DiseaseLocus> ModelLociArray;

}
}
#endif //DISEASELOCI

