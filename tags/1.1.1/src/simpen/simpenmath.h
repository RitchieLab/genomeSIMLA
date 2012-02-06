/*

  BcwMath.h - Bill White - 11/02/01

  Bill's general math routines. See BcwMath.C file for function descriptions.
*/

#ifndef _EVOPEN_MATH_H_
#define _EVOPEN_MATH_H_



#include <map>
#include <iostream>
#include<sstream>

#include "types.h"
#include "tree.hh"
#include "node.h"

namespace SimPen {



using namespace std;


#define ZERO_THRESHOLD 1e-17


void CalculateMarginals(double * tableValue2, double * marginals2, vector< vector <double> > & geno_freqs,
  vector< tree<node *> > & locus_trees);

void construct_trees(vector< tree<node *> > & locus_trees,
  unsigned int num_loci);
void nested_tree_loops(int max_depth, unsigned int * lower, unsigned int *upper,
  unsigned int * locus_order, unsigned int * mult_indexes,
  tree<node *> & tr);
int get_index_value(unsigned int * mults, unsigned int * indexes,
  int cur_depth);

void createTrees(vector< vector< tree<node *> > >& subModelTrees, int numLoci,
  vector< tree<node *> > & existingTrees);
 
void checkSubModels(double* penValues, int numLoci, vector<int> lociList,
  map<string, resultInfo> & resultMap, vector< vector <double> > & geno_freqs,
  vector< vector< tree<node *> > >& subModelTrees);
  
double* collapsePenTable(double* penValues, vector< vector <double> > & geno_freqs,
  vector<int> & lociList, int excludedLocIndex);
  
string itos(int number);

double calculate_tree(tree<node *>::iterator child, tree<node *> & tr,
   vector< vector<double> > & geno_freqs,double * table_values);
double calculate_prob(tree<node *>::iterator location, tree<node *> & tr,
   vector< vector<double> > & geno_freqs,double * table_values);
void clear_tree(tree<node *> & tr);

double CalculateVariance(double* values, int numValues);

double CalculateHeritability(tree<node *> & herit_tree,
          double * marginalProbability, int num_loci, int num_genos_per_locus,
          vector< vector <double> > & geno_freqs, double * table_values, double * marginal_mean);
double calculate_herit(tree<node *>::iterator child, tree<node *> & herit_tree,
  vector< vector <double> > & geno_freqs, double marginal_mean, double * table_values);
double calculate_marginal_mean(vector< vector <double> > & geno_freqs, double * marginalProbs,
   int num_loci, int num_genos_per_locus);


}
#endif
