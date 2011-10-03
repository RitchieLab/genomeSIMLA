/*
  EvolvePenetrance.h - Bill White - 11/01/01

  Header file for GA-evolved Penetrance Tables Project
*/

#ifndef _SIMPEN_H_
#define _SIMPEN_H_

#include "types.h"

#include "tree.hh"
#include "node.h"

#include "simpreferences.h"

#include <ga/ga.h>
#include <ga/GARealGenome.h>

// include math.h for FLT_MAX on Debian
#include <math.h>
#include "simpenmdr.h"

// maximums
#define MAX_STRING 2048
#define MAX_LOCI 1024

// ERRORS
#define COMMAND_LINE_ERROR -1
#define CONFIG_FILE_ERROR -2
#define ANT_COLONY_INIT_ERROR -3
#define DATA_FILE_ERROR -4
#define CV_INIT_ERROR -5
#define ARRAY_BOUNDS_ERROR -6
#define DEBUG_ERROR -98
#define FATAL_ERROR -99

// booleans
#define FALSE 0
#define TRUE 1

// epsilon is the maximum misclassification error
// to consider a "hit"
#define EPSILON 0.000001

namespace SimPen {



//int run_simpen(int argc, char** argv);
float run_simpen(const char *filename, const char *output, Simulation::StatusModel::ModelLociArray &loci, int seed, bool override);

// function prototypes
float Objective(GAGenome& g);
void AddModelToBest(ofstream & best_model_file, GARealGenome& genome, int num_loci,
     const vector<string> & labels);
void CreateLabels(int num_loci, vector<string> & labels);
void header_nested_loops( unsigned int * lower, unsigned int *upper, int max_depth,
     vector<string> & labels, vector< vector<string> > & genotypes);
void summarize_scores(GARealGenome& g, ofstream & best_model_file, SimPreferences & pref,
  vector< vector< tree<node *> > >  & subModelTrees);
void header(string progname, string versiondate);
void SubModelCheck(double * tableValue, unsigned int genome_length,  ofstream & gmFile, SimPreferences & pref,
  vector< vector< tree<node *> > >  & subModelTrees);
  
void show_score_update(GARealGenome& g);
void calc_score(GARealGenome& genome, double * marginalProbability, double & marginAvg,
  double & heritability, double & tableVariance, double & marginVariance, double & oddsRatio,
  bool oddsCalc, bool heritCalc);
double calc_fitness(double tableVariance, double marginAvg, double marginVariance,
  double heritability, double oddsRatio);
  
string CreateModelFile(string base_file, int num_loci, GARealGenome& genome,
  const vector<string> & labels, int modelNum, Simulation::StatusModel::ModelLociArray &loci);

const string program_name = "simpen";
const string version_date = "11/17/05";







}

#endif
