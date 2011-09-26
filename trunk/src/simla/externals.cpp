#include "externals.h"

//*********************************

void generate_genotypes(int loci_ex, double *freq, int n, int *list)
{
	int i,j;


	double a1, a2;
	int total;
	int type;
	int multiplier;
	for(i = 0; i < n; i++){
		total = 0;
		multiplier = (int)powf(3,(float)(loci_ex - 1));					// 3 distinct genotype
		for(j = 0; j <= (loci_ex - 1); j++){
			a1 = r.ran1();
			a2 = r.ran1();
			if((a1 < freq[j]) && (a2 < freq[j])){
				type = 2;				// (1,1) genotype
			}
			else if((a1 >= freq[j]) && (a2 >= freq[j])){
				type = 0;				// (2,2) genotype
			}
			else{
				type = 1;				// (1,2) genotype
			}
			total += (multiplier * type);
			multiplier /= 3;
		}
	 	list[i] = total;
	}
}

//*********************************

//*********************************

//*********************************

//*********************************

//*********************************

//*********************************
