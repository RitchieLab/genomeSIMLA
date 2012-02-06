#ifndef SIMPEN_TYPES
#define SIMPEN_TYPES


#include <vector>


namespace SimPen {

using namespace std;

// min and max values for real genes
const double MIN_GENE_VALUE = 0.0;
const double MAX_GENE_VALUE = 1.001;
const double INCREMENT = 0.001;

// # of genotypes for each locus
const int NUM_GENOS_PER_LOCUS = 3;

const double STOPPING_VALUE = 1.0;

const double VARIANCE_THRESHOLD = 1e-5;

struct resultInfo{
  vector<int> lociList;
  double * margPenetrances;
  double marginVariance;
  double tableVariance;
  double herit;
  double marginAvg;
  double oddsRatio; 
};

}

#endif
