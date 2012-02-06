// simpenMDR.C

#include "simpenmdr.h"
#include <iostream>

namespace SimPen {

void genotype_freqs_nested_loops( unsigned int * lower, unsigned int *upper, 
  int max_depth, double * geno_freq_list, vector< vector <double> > & geno_freqs);

// Functions for creating the MDR model and calculating
// the odds ratio
// returns an array that lists the genotype frequencies 
// corresponding to each position in the array
// Args: number of loci
//       number of genotypes per locus
//       2-D vector of genotype frequencies
//       array long enough to hold all the frequencies
// Ret:  array listing frequencies that will match 
//       the genome order for penetrance values
double * set_genotype_freqs(int num_loci, int num_genos,
  vector< vector <double> > & geno_freqs, double * geno_freq_list){
  
  unsigned int * lower_index = new unsigned int [num_loci];
  unsigned int * upper_index = new unsigned int [num_loci];
  
  for(int i=0; i<num_loci; i++){
    lower_index[i] = 0;
    upper_index[i] = num_genos;
  }
   
  genotype_freqs_nested_loops(lower_index, upper_index, num_loci, 
    geno_freq_list, geno_freqs);  
  delete [] lower_index;
  delete [] upper_index;
  
  return geno_freq_list;
}

// performs nested loops for setting genotype frequencies for each 
// combination in the model
void genotype_freqs_nested_loops( unsigned int * lower, unsigned int *upper, 
  int max_depth, double * geno_freq_list, vector< vector <double> > & geno_freqs){

   unsigned int * indexes = new unsigned int[max_depth];
   int cur_depth = 0;
   int geno_freq_index = 0;
   indexes[cur_depth] = lower[cur_depth];

   while(1){
      if( indexes[cur_depth] < upper[cur_depth] ) {
         if( cur_depth == max_depth - 1 ) {
           // calculate frequency here
           geno_freq_list[geno_freq_index] = 1;
           for(int i =0; i<=cur_depth; i++){
              geno_freq_list[geno_freq_index] *= geno_freqs[i][indexes[i]];
            }
           geno_freq_index++;
           indexes[cur_depth]++;
         }
        else{
          cur_depth++;
          indexes[cur_depth]=lower[cur_depth];
        }
      }
      else{
        if(cur_depth > 0){
          indexes[--cur_depth]++;
        }
        else{
          break;
        }
      }
   }
   delete [] indexes;
}



// Calculates and returns the odds ratio based on the penetrance
// table 
// Arg:  array of penetrance values
//       array of genotype frequencies that match order of penetrance values
//       number of cells in the penetrance table
// Ret: odds ratio
double calculate_odds_ratio(double * penetranceValues, double * geno_freq_list,
  int num_cells){
  
  double unaffectedHigh=0.0, unaffectedLow=0.0, affectedHigh=0.0, affectedLow=0.0;
  double * unaffectedFreqs = new double [num_cells];
  double * affectedFreqs = new double [num_cells];
  double unaffectedSum=0.0, affectedSum=0.0;
  
  // set cells in model as high 1 or low 0
  for(int cellIndex=0; cellIndex < num_cells; cellIndex++){
    unaffectedFreqs[cellIndex] = (1-penetranceValues[cellIndex]) * geno_freq_list[cellIndex];
    affectedFreqs[cellIndex] = penetranceValues[cellIndex] * geno_freq_list[cellIndex];
    unaffectedSum += unaffectedFreqs[cellIndex];
    affectedSum += affectedFreqs[cellIndex];
  }
  
  // determine the fraction of individuals that will
  // have each genotype for affected and unaffected
  for(int cellIndex=0; cellIndex < num_cells; cellIndex++){
    unaffectedFreqs[cellIndex] = unaffectedFreqs[cellIndex] / unaffectedSum;
    affectedFreqs[cellIndex] = affectedFreqs[cellIndex] / affectedSum;
    // This is a high incidence cell
    if(affectedFreqs[cellIndex] > unaffectedFreqs[cellIndex]){
      unaffectedHigh += unaffectedFreqs[cellIndex];
      affectedHigh += affectedFreqs[cellIndex];
    }
    else{
      unaffectedLow += unaffectedFreqs[cellIndex];
      affectedLow += affectedFreqs[cellIndex];
    }
  }
  
  double oddsRatio;
  // set oddsRatio to very high number when it would be infinite
  if(unaffectedHigh == 0 || affectedLow ==0){
    oddsRatio = 100000;
  }
  else{
    oddsRatio = affectedHigh * unaffectedLow / (unaffectedHigh * affectedLow);
  }
  
  delete [] unaffectedFreqs;
  delete [] affectedFreqs;

  return oddsRatio;
}

}
