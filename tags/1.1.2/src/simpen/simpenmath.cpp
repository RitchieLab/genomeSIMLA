/*
simpenMath.C

Scott Dudek 7/19/04

Used in generating models.

Adapted from Bill White.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "simpenmath.h"
#include "simpenmdr.h"

namespace SimPen {


/**************************************************************************
 * CalculateVariance - Calculate the variance of a vector of doubles
 *
 * Input:  vector of doubles, number of values in vector
 * Output: variance
 **************************************************************************/
double CalculateVariance(double* values, int numValues)
{
        double sum = 0.0;
        int i;

        // calculate average
        for(i=0; i < numValues; i++) {
                sum += values[i];
        }
        double average = sum / numValues;

        // calculate and return variance
        double SSE = 0.0;
        for(i=0; i < numValues; i++) {
                SSE += ((values[i] - average) * (values[i] - average));
        }
        return (SSE / (numValues - 1));
}


/**************************************************************************
 * CalculateMarginals - Calculate the marginal penetrance values for a

 *                      frequencies
 *
 * Input:  2-locus model as 3x3 matrix, pointer to vector of doubles to
 *         return marginal values, allele frequencies
 * Output: vector of marginal values
 **************************************************************************/
void CalculateMarginals(double pValues[],
   double marginalProbability[], vector< vector <double> > & geno_freqs,
   vector< tree<node *> > & locus_trees)
{

   // now loop through and calculate all the other trees
   // to determine the marginal probabilities
   tree<node *>::iterator root, child;

   for(unsigned int i=0; i<locus_trees.size(); i++){
      //Assign table values to leaves of tree for calculating
      //marginal probability
      root = locus_trees[i].begin();

      // Calculate the 3 marginal probabilities for each tree
      for(int j=0; j<3; j++){
        child = locus_trees[i].child(root, j);
        marginalProbability[i*3+j] = calculate_tree(child,locus_trees[i], geno_freqs, pValues);
      }
   }
}


/**************************************************************************
* construct_trees:  constructs trees for use in calculations
*
*
**************************************************************************/
void construct_trees(vector< tree<node *> > & locus_trees,
   unsigned int num_loci){
   unsigned int * mult_indexes = new unsigned int[num_loci];
   unsigned int * locus_order = new unsigned int[num_loci];
   unsigned int * lower = new unsigned int[num_loci];
   unsigned int * upper = new unsigned int[num_loci];

   for(unsigned int i = 0; i<num_loci; i++){
      lower[i] = 0;
      upper[i] = 3;
      mult_indexes[i] = (unsigned int)(pow(3,double(num_loci-i-1)));
      locus_order[i] = i;
   }

   int last_index = num_loci-1;

    for(unsigned int i=0; i<num_loci; i++){
      tree<node *> tr;
      // after first run move multipliers and locus index
      // to track which ones are currently being calculated
      // Starts with A then goes to B, C etc.
     if(i > 0){
        unsigned int temp_l = locus_order[0];
        unsigned int temp_m = mult_indexes[0];
        for(int k=0; k<last_index; k++){
          locus_order[k] = locus_order[k+1];
          mult_indexes[k] = mult_indexes[k+1];
        }

        locus_order[last_index] = temp_l;
        mult_indexes[last_index] = temp_m;
      }
      // pass tree into loops along with indexes and action
      nested_tree_loops(num_loci, lower, upper, locus_order,
        mult_indexes, tr);
      locus_trees.push_back(tr);
   }

   delete [] upper;
   delete [] lower;
   delete [] locus_order;
   delete [] mult_indexes;
}


/**************************************************************************
* nested_tree_loops - returns tree for evaluating penetrances
* Input - current locus index,
*         number of loci,
*         array of penetrance table values
* Ret   - tree
 **************************************************************************/
void nested_tree_loops(int max_depth, unsigned int * lower, unsigned int *upper,
   unsigned int * locus_order, unsigned int * mult_indexes, tree<node *> & tr){


   unsigned int * indexes = new unsigned int[max_depth];
   int cur_depth = 0;

   tree<node * >::iterator current_node, top;

   top=tr.begin();
   current_node=tr.insert(top, new node());

   indexes[cur_depth] = lower[cur_depth];

   while(1){
      if( indexes[cur_depth] < upper[cur_depth] ) {
         if( cur_depth == max_depth - 1 ) {
           tr.append_child(current_node, new node(locus_order[cur_depth],
            indexes[cur_depth], get_index_value(mult_indexes, indexes, cur_depth)));
           indexes[cur_depth]++;
         }
        else{
          current_node = tr.append_child(current_node, new node(locus_order[cur_depth],
            indexes[cur_depth]));
          cur_depth++;
          indexes[cur_depth]=lower[cur_depth];
        }
      }
      else{
        if(cur_depth > 0){
          current_node = tr.parent(current_node);
          indexes[--cur_depth]++;
        }
        else{
          break;
        }
      }
   }
   delete [] indexes;
 }


 /**************************************************************************
* get_index_value - returns index for that node in
* Input - integer array of multiplier for determining location in array,
*         number of loci,
*         integer array of indexes (representing the theoretical multi-
*         dimensional array being scanned,
* Ret   - index of one-dimensional array for the table value
 **************************************************************************/
 int get_index_value(unsigned int * mults, unsigned int * indexes,
      int cur_depth){
   int table_index = 0;

   for(int i =0; i<=cur_depth; i++){
      table_index += indexes[i] * mults[i];
   }
   return table_index;
 }


/**************************************************************************
* clear_tree
*
*
**************************************************************************/
void clear_tree(tree<node *> & tr){
      tree<node *>::iterator iter;
      for(iter=tr.begin(); iter!=tr.end(); iter++){
        delete (*iter);
      }
      tr.clear();
}


/**************************************************************************
* calculate_tree
*
*
**************************************************************************/
double calculate_tree(tree<node *>::iterator start,tree<node *> & tr,
   vector< vector<double> > & geno_freqs, double * table_values){

   double marginalProb = 0.0;
   tree<node *>::iterator child;
   unsigned int num_children = tr.number_of_children(start);

   for(unsigned int i=0; i<num_children; i++){
      child = tr.child(start, i);
      marginalProb += calculate_prob(child, tr, geno_freqs, table_values);
   }


   return marginalProb;
}

 /**************************************************************************
* calculate_prob - recursive algorithm for determining portion of marginal
*   probability of sub tree
*
*
*
**************************************************************************/
double calculate_prob(tree<node *>::iterator location, tree<node *> & tr,
   vector< vector<double> > & geno_freqs, double * table_values){

   // when no children below return value stored in leaf multiplied by
   // the frequency
   if(tr.number_of_children(location) == 0){
      return  table_values[(*location)->index] * geno_freqs[(*location)->locus][(*location)->genotype_num];
   }

   double total = 0.0;
   unsigned int num_children = tr.number_of_children(location);
   tree<node *>::iterator child;
   for(unsigned int i=0; i<num_children; i++){
      child = tr.child(location, i);
      total += calculate_prob(child, tr, geno_freqs, table_values);
   }
   return total * geno_freqs[(*location)->locus][(*location)->genotype_num];
}


/**************************************************************************
* CalculateHeritability - returns heritability
* Input - tree constructed with each level being one genotype,
*         double array of marginal Probabilities,
*         number of loci
*         number of genotypes per locus
*         2-d vector of genotype frequencies
* Ret   - heritability
 **************************************************************************/
double CalculateHeritability(tree<node *> & herit_tree,
          double * marginalProbability, int num_loci, int num_genos_per_locus,
          vector< vector <double> > & geno_freqs, double * table_values, double * marginal_mean){

   *marginal_mean = calculate_marginal_mean(geno_freqs, marginalProbability,
      num_loci, num_genos_per_locus);

   tree<node *>::iterator root, child;
   // iterate through tree to calculate the heritability
    root = herit_tree.begin();
    double VI = 0.0;
   for(int j=0; j<num_genos_per_locus; j++){
        child = herit_tree.child(root, j);
        VI += calculate_herit(child, herit_tree, geno_freqs, *marginal_mean, table_values);
   }

   return (VI/(*marginal_mean * (1.0 - *marginal_mean)));
}


/**************************************************************************
* calculate_herit - recursively traverses tree and calculates
*                         numerator for heritability formula
* Input - iterator pointing to node to check
*         tree constructed with each level being one genotype,
*         2-d vector of genotype frequencies
*         marginal mean
* Ret   - VI in heritability calculation
 **************************************************************************/
double calculate_herit(tree<node *>::iterator location, tree<node *> & herit_tree,
   vector< vector <double> > & geno_freqs, double marginal_mean, double * table_values){
    // when no children below return value stored in leaf multiplied by
   // the frequency
   if(herit_tree.number_of_children(location) == 0){
      return (table_values[(*location)->index] - marginal_mean) *
        (table_values[(*location)->index] - marginal_mean) *
        geno_freqs[(*location)->locus][(*location)->genotype_num];
   }

   double total = 0.0;
   unsigned int num_children = herit_tree.number_of_children(location);
   tree<node *>::iterator child;
   for(unsigned int i=0; i<num_children; i++){
      child = herit_tree.child(location, i);
      total += calculate_herit(child, herit_tree, geno_freqs, marginal_mean, table_values);
   }
   return total * geno_freqs[(*location)->locus][(*location)->genotype_num];
}


/**************************************************************************
* calculate_marginal_mean - returns penetrance value for indicated location
* Input - 2-d vector of genotype frequencies
*         double array of marginal probabilities
*         number of loci
*         number of genotypes per locus
* Ret   - mean of marginal probabilities
 **************************************************************************/
double calculate_marginal_mean(vector<vector <double> > & geno_freqs, double * marginalProbs,
   int num_loci, int num_genos_per_locus){
   double marginSum = 0.0;

   int i,j, marginProbCount=0;;
   for(i=0; i<num_loci; i++)
      for(j=0; j<num_genos_per_locus; j++)
         marginSum += geno_freqs[i][j] * marginalProbs[marginProbCount++];

   return marginSum / num_loci;
}

/***************************************************************************
* createTrees - Creates 2-D vector where first index is number of Loci
*               in the model.  These trees are used in determining stats
*               for sub models
*
****************************************************************************/
void createTrees(vector< vector< tree<node *> > >& subModelTrees, int numLoci,
  vector< tree<node *> > & existingTrees){
  
  vector< tree<node *> > emptyVector;
  for(int lociCount=0; lociCount < numLoci; lociCount++){
    subModelTrees.push_back(emptyVector);
  }
  
  // put trees into vector for existing number of loci
  subModelTrees.push_back(existingTrees);

  // now create trees for each of the lesser sets
  // down to 2
  for(int lociCount=numLoci-1; lociCount >= 2; lociCount--){
    construct_trees(subModelTrees[lociCount], lociCount);
  }
}

/***************************************************************************
* checkSubModels - checks sub models of the list that is passed in
*                - works recursively to track the results
*
****************************************************************************/
void checkSubModels(double* penValues, int numLoci, vector<int> lociList,
  map<string, resultInfo> & resultMap, vector< vector <double> > & geno_freqs,
  vector< vector< tree<node *> > > & subModelTrees){

  resultInfo newResult;
  int newNumLoci = numLoci-1;

  // drop each locus in list and check remainder for model
  for(int currLocusIndex = 0; currLocusIndex < numLoci; currLocusIndex++){
    vector<int> newLociList;
    string combinationStr;
    // delete entry of current locus
    for(int counter=0; counter < numLoci; counter++){
      if(counter != currLocusIndex){
        combinationStr += itos(lociList[counter]) + ".";
        newLociList.push_back(lociList[counter]);
      }
    }
    // check to see if the combination matches any existing results
    // currLocusIndex indicates which index should be dropped
    if(resultMap.find(combinationStr) == resultMap.end()){
      // when combination not found proceed to create new values and test
      double * newValues = collapsePenTable(penValues, geno_freqs, lociList, currLocusIndex);
      // have to put correct genotype frequencies in 2-D array for calculating
      // marginal penetrances
      vector< vector <double> > newGenoFreqs;
      for(int locus=0; locus<newNumLoci; locus++){
        newGenoFreqs.push_back(geno_freqs[newLociList[locus]]);
      }
      
      // create new genotype frequencies for use in Odds Ratio calculation
      unsigned int genotype_freq_length = 
        (unsigned int)(pow(float(NUM_GENOS_PER_LOCUS), newNumLoci));
      double * new_genotype_freqs = new double[genotype_freq_length];   
      new_genotype_freqs = set_genotype_freqs(newNumLoci, NUM_GENOS_PER_LOCUS,
        geno_freqs, new_genotype_freqs);
      
      newResult.lociList = newLociList;
      newResult.margPenetrances = new double[newNumLoci*NUM_GENOS_PER_LOCUS];
      CalculateMarginals(newValues, newResult.margPenetrances, newGenoFreqs, 
        subModelTrees[newNumLoci]);
      
      newResult.marginVariance = CalculateVariance(newResult.margPenetrances, 
        newNumLoci*NUM_GENOS_PER_LOCUS);
        
      newResult.tableVariance = CalculateVariance(newValues, int(pow(3,double(newNumLoci))));

      newResult.herit = CalculateHeritability(subModelTrees[newNumLoci][0],
        newResult.margPenetrances, newNumLoci, NUM_GENOS_PER_LOCUS, newGenoFreqs, newValues, &newResult.marginAvg);
      
      newResult.oddsRatio = calculate_odds_ratio(newValues,new_genotype_freqs, 
        genotype_freq_length);

      resultMap[combinationStr] = newResult;   
      
      // make recursive call using passing list of loci and new table
      if(newNumLoci > 2){
        checkSubModels(newValues, newNumLoci, newLociList, resultMap,
          geno_freqs, subModelTrees);
      }
      
      // use new penetrance values along with trees for that size 
      // of table to calculate marginal Penetrances      
      delete [] newValues;
      delete [] new_genotype_freqs;
    }
  }
}

/***************************************************************************
* itos - Convert integer into string
*
****************************************************************************/
string itos(int number){
  stringstream oss;
  oss << number;
  return oss.str();
}


/***************************************************************************
* collapsePenTable - reduces penetrance table by removing a locus
*                    81 cell 4-locus table becomes a 27 cell 3-locus table
*
****************************************************************************/
double* collapsePenTable(double * penValues, vector< vector <double> > & geno_freqs,
  vector<int> & lociList, int excludedLocIndex){
  
  unsigned int numLoci = lociList.size();
  
  double * newPenValues = new double[int(pow(3,double(numLoci-1)))];

  // these are used to multiply indexes to determine which values of
  // the original penetrance table to use

  unsigned int * lower = new unsigned int[numLoci];
  unsigned int * upper = new unsigned int[numLoci];  
  unsigned int * locus_order = new unsigned int[numLoci];
  unsigned int * mult_indexes = new unsigned int[numLoci];
  unsigned int * indexes = new unsigned int[numLoci];
  
  int lastIndex = numLoci-1;
  int max_depth = numLoci;
  
  // initialize the arrays
  for(unsigned int i=0; i<numLoci; i++){
    mult_indexes[i] = (unsigned int)(pow(3,double(numLoci-i-1)));
    lower[i] = 0;
    upper[i] = 3;
    locus_order[i] = lociList[i];
  }
  
  // now swap array positions if excludedLocus is not the last one
  if(excludedLocIndex != lastIndex){
    int tempLocOrder = locus_order[excludedLocIndex];
    int tempMultIndex = mult_indexes[excludedLocIndex];
    // place excluded locusIndex in last position and shift rest left
    for(int i=excludedLocIndex; i<lastIndex; i++){
      locus_order[i] = locus_order[i+1];
      mult_indexes[i] = mult_indexes[i+1];
    }
    locus_order[lastIndex] = tempLocOrder;
    mult_indexes[lastIndex] = tempMultIndex;
  }

  int cur_depth = 0;
  indexes[cur_depth] = lower[cur_depth];
  int newPenTableIndex=0;
  double tableCell = 0.0;
  // loop goes through all positions and consolidates the values
  // into a new list that has one less locus
  while(1){

    if( indexes[cur_depth] < upper[cur_depth] ) {
      if( cur_depth == max_depth - 1 ) {
        // at lowest level multiply variable by frequency
        int oldPenIndex = get_index_value(mult_indexes, indexes, cur_depth);
        tableCell +=  penValues[oldPenIndex] * geno_freqs[locus_order[lastIndex]][indexes[cur_depth]];
        if(++indexes[cur_depth] == upper[cur_depth]){
          newPenValues[newPenTableIndex++] = tableCell;
          tableCell = 0.0;
        }
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

  delete [] upper;
  delete [] lower;
  delete [] locus_order;
  delete [] mult_indexes;
  delete [] indexes;

  return newPenValues;
}


}

 
