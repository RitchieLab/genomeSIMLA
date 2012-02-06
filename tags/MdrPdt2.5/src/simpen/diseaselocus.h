#ifndef DISEASELOCI
#define DISEASELOCI

#include <math.h>
#include <vector>

namespace SimPen {

using namespace std;


struct DiseaseLocus {
	int chromosome;					///<The chromosome the locus is located on
	int locusIdx;					///<The index of the locus
	float alFreq1;					///<Allele frequency
	float alFreq2;					///<Allele frequency
	
	DiseaseLocus() : chromosome(0), locusIdx(0), alFreq1(0.0), alFreq2(0.0) { }
	DiseaseLocus(int c, int i, float f1, float f2) : chromosome(c), locusIdx(i), alFreq1(f1), alFreq2(f2) { }

	/**
	 * @brief Returns that the frequencies are valid (sum to 1.0)
	 */
	bool Verify() {  return fabs(1.0 - (alFreq1 + alFreq2)) < 0.001; }

	void Report() { cout<<"\t"<<chromosome<<":"<<locusIdx<<" "<<alFreq1<<" "<<alFreq2<<" "<<fabs(1.0 - (alFreq1 + alFreq2))<<"\n"; }


};


typedef std::vector<DiseaseLocus> LocusArray;

}

#endif //DISEASELOCI

