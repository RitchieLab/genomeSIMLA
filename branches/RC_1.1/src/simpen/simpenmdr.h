// simpenMDR.h

#ifndef _SIMPENMDR_H_
#define _SIMPENMDR_H_

#include <vector>

namespace SimPen {

using namespace std;

double * set_genotype_freqs(int num_loci, int num_genos,
  vector< vector <double> > & geno_freqs, double * geno_freq_list);
  
double calculate_odds_ratio(double * penetranceValues, double * geno_freq_list,
  int num_cells);

}

#endif
