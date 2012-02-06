/*

   11/15/04 - smd
   Added check for models involving loci that make up overall model
   

   Scott Dudek (smd)- 7/1/04
   Adapted from Bill White and Marylyn Ritchie's code for evolving and
   evaluating penetrances models.

                EvolvePenetrance.C - Bill White - 10/29/01

                adapted from:
                ga_models.C - Bill White - 1/7/99

                adapted from:
                ex23.C
                mbwall 5jan96
                Copyright (c) 1995-1996  Massachusetts Institute of Technology
*/

#include <fstream>
#include <iomanip>
//#include <unistd.h>
#include <sstream>

#include <iostream>
#include <math.h>

#include "simpen/prefexcept.h"
#include "simpen/simpen.h"
#include "simpen/simpenmath.h"


// GAlib objects       
#include <ga/ga.h>
#include <ga/GARealGenome.h>
#include <ga/GARealGenome.C>			//o.O - This is required due to a trick MIT is doing
										//to get the alleles automatically instanced (and gcc doesn't
										//do this. I think we can use: __attribute__((weak))__ to 
										//fix this, but it will require changing the code in galib. 
										//This goes on the distant todo list


#include "utility/utility.h"

namespace SimPen {

using namespace Utility;
using namespace Simulation::StatusModel;


vector< vector <double> > geno_freqs;
vector< tree<node *> > locus_trees;
unsigned int lociCount			= 0;
double MARGINS_VAR_THRESHOLD 	= 0.0;
double TABLE_VAR_THRESHOLD 		= 0.0;
double HERITABILITY_THRESHOLD 	= 0.0;
double TARGET_PEN				= 0.0;
double TARGET_PEN_DIFF 			= 0.003;
double ODDS_TARGET				= 0.0;

bool CALCODDS 					= false;
bool CALCHERIT 					= false;
double * GENO_FREQUENCY_LIST 	= NULL;

int MARGIN_WEIGHT				= 0;
int HERIT_WEIGHT				= 0;
double ODDS_WEIGHT				= 0.0;



// -------------------------------- M A I N -----------------------------------
float run_simpen(const char *filename, const char *output, ModelLociArray &loci, int seed, bool override) {
	float bestFitness =0.0;

	//using namespace SimPen;
	using namespace std;

	string controlfile, output_file = output, pen_file, datasimFileName;
  
	pen_file = output_file;
	output_file += ".out";

	SimPreferences pref;

	cout<<"Attempting to read model configuration file: "<<filename<<"\n";
	try{
		pref.read_config(filename);
	}
	catch(PrefExcept & pe){
		cerr << pe << endl << endl;
		return 0.0;
		//throw FileIOError(filename, 0);
	}


	lociCount = loci.size();
	//OK, here we will set up the disease loci according to the parameters (over-ruling whatever
	//might have been found in the config file.
	for (uint i=0; override && i<lociCount; i++) {
		pref.AddDiseaseLocus(loci[i]);
		//pref.AppendFrequencies(loci[i].alFreq1, loci[i].alFreq2);
	}


	// See if we've been given a seed to use (for testing purposes).  When you
	// specify a random seed, the evolution will be exactly the same each time
	// you use that seed number.
	int update_interval = pref.get_update_interval();

	if (override && seed > 0) 
		pref.SetSeed(seed);

	GARandomSeed(seed);
	
	// -------------------------- S E T U P  G A --------------------------------
	
	GARealAlleleSet alleles(MIN_GENE_VALUE, MAX_GENE_VALUE, INCREMENT);
	// each chomosome has  # of genotypes per loci ^ number of loci
	// for 2 snps total is 3 ^ 2 = 9
	unsigned int total_length = (unsigned int)(pow(float(NUM_GENOS_PER_LOCUS), (int)lociCount));
	

	if (override)
		pref.ResetModel();
	//OK, here we will set up the disease loci according to the parameters (over-ruling whatever
	//might have been found in the config file.
	for (uint i=0; override && i<lociCount; i++) {
		pref.AddDiseaseLocus(loci[i]);
		//pref.AppendFrequencies(loci[i].alFreq1, loci[i].alFreq2);
	}

	//pref.set_loci(lociCount);
	//vector< vector<double> > all_freqs = pref.get_freqs();


	// set up the allele frequencies from the SimPreferences file
	uint i;
	geno_freqs.clear();
	cout<<"\n\nGenotype Frequencies: \n";

	cout<<"\t"<<setw(5)<<"loc"<<" "<<setw(22)<<"Allele Freq.      "<<setw(25)<<"Genotype Freq.\n";

	for(i=0; i<lociCount; i++){
		//Hm, this looks like it needs the genotypes to go in backwards
		DiseaseLocus &l = pref.GetDiseaseLocus(i);
		geno_freqs.push_back(vector<double>(3,0));
		geno_freqs[i][0] = l.alFreq1 * l.alFreq1;
		geno_freqs[i][1] = l.alFreq1 * l.alFreq2 * 2;
		geno_freqs[i][2] = l.alFreq2 * l.alFreq2;

		cout<<"\t"<<setw(5)<<i<<" "<<setw(10)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<l.alFreq1<<" ";
		cout<<setw(10)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<l.alFreq2<<" : ";
		cout<<setw(10)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<geno_freqs[i][0]<<" ";
		cout<<setw(10)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<geno_freqs[i][1]<<" ";
		cout<<setw(10)<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<geno_freqs[i][2]<<"\n";
 	}

	if (GENO_FREQUENCY_LIST)
		delete[] GENO_FREQUENCY_LIST;
  	GENO_FREQUENCY_LIST = new double[total_length];
  	set_genotype_freqs(lociCount, NUM_GENOS_PER_LOCUS, geno_freqs, GENO_FREQUENCY_LIST);

  

	MARGINS_VAR_THRESHOLD = pref.get_marg_var();
	
	MARGIN_WEIGHT = pref.get_marg_weight();
	HERIT_WEIGHT = pref.get_herit_weight();
	ODDS_WEIGHT = pref.get_odds_weight();
	ODDS_TARGET = pref.get_odds_ratio();
  
	if(ODDS_TARGET > 0){
		CALCODDS = true;
	}
	if(pref.get_tablevar_on())
    	TABLE_VAR_THRESHOLD = pref.get_table_var();
	else
		TABLE_VAR_THRESHOLD = 0.0;
	HERITABILITY_THRESHOLD = pref.get_herit();	
	if(HERITABILITY_THRESHOLD > 0){
    	CALCHERIT = true;
	}
	TARGET_PEN = pref.get_target_penetrance();
	
	GARealGenome genome(total_length, alleles, Objective);
	GAUniformSelector selector;
	
	//GASimpleGA ga(genome);
	GADemeGA ga(genome);
	ga.maximize();                // we want to maximize the objective function
	
	ga.scaling(GANoScaling());                // set the scaling method to our sharing
	ga.selector(selector); // set selection method
	ga.nPopulations(pref.get_demes());
	ga.populationSize(pref.get_pop_size());        // how many individuals in the population
	ga.nGenerations(pref.get_gen());                // number of generations to evolve
	ga.pMutation(pref.get_mutate());                // likelihood of mutating new offspring
	ga.pCrossover(pref.get_cross());                // likelihood of crossing over parents
	
	ga.selectScores(GAStatistics::AllScores);



	construct_trees(locus_trees, lociCount);  // constructs trees for use in calculations
	vector< vector< tree<node *> > >  subModelTrees;
	createTrees(subModelTrees, lociCount, locus_trees);
	
	// evolve MAX_MODEL_LOOPS models
	ofstream gmFile(output_file.c_str(), ios::out);
	vector<string> labels;
	CreateLabels(lociCount, labels);
	for(int curModelLoop=0; curModelLoop < pref.get_num_models(); curModelLoop++) {
	
		// --------------------------- R U N  G A ---------------------------------
		
		printf("Model: %5d of %5d ...Starting evolution...\n\n",
			curModelLoop+1, pref.get_num_models());
		if(pref.get_verbose()){
			cout << setw(10) << "Count" << setw(14) << setfill(' ') << "MarginalVar";
			cout << setw(14) << "TableVar" << setw(14) << "Heritability";
			cout << setw(12) << "MarginAvg" << setw(12) << "Odds Ratio" << setw(10) << "Fitness" << endl;
			cout << "---------------------------------------------------------------------------------------------" << endl;
		}
		else{
			cout << setw(10) << "Count" << setw(10) << "Fitness" << endl;
			cout << "---------------------" << endl;
		}
		int printCounter = 0;
		ga.initialize();
		int converged = 0;
		while(!ga.done() && !converged) {
			// go to next generation
			ga.step();
		
			// print best of generation
			if( printCounter % update_interval == 0) {
				stringstream ss;
				ss << "[" << printCounter << "]";
				cout << setw(10) << ss.str();
				if(pref.get_verbose()){
					GARealGenome& g = (GARealGenome&)ga.statistics().bestIndividual();
					show_score_update(g);
				}
				cout << setw(10) <<  setprecision(6) << ga.statistics().bestIndividual().score() << endl;             
			}
		
			if(ga.statistics().bestIndividual().score() >= STOPPING_VALUE) {
				converged = 1;
			}
			else {
				printCounter++;
			}
		}
	
		// ------------------- P R O C E S S   G A   B E S T ----------------------
	
		// save this model to the PenetranceTable
		GARealGenome& genome = (GARealGenome&) ga.statistics().bestIndividual();
		if(genome.score() <= 1.0) {
	
			cout << endl << "Model: " << setw(5)  << curModelLoop + 1 << " of " <<
				setw(5) << pref.get_num_models() << " ...Model with best fit being saved...\n\n";
				gmFile << "-----------------------------------" << endl;    
				gmFile << "Model #" << curModelLoop+1 << endl;
				gmFile << "-----------------------------------" << endl;
			
			AddModelToBest(gmFile, genome, lociCount, labels);
			if(pref.get_verbose() || pref.get_submodels_on())
				summarize_scores(genome, gmFile, pref, subModelTrees);
		
		
			// for each model create sample datasimulation file
			CreateModelFile(pen_file, lociCount, genome, labels, curModelLoop+1, loci);
			
			//Eventually, we need to set this according to the true winner
			bestFitness = genome.score();
		}
		else {
			cout << "[" << setw(5) << curModelLoop + 1<< "/" <<
				setw(5) << pref.get_num_models() << "] Best-of-run score["
				<< genome.score() << "] above threshold [1].\n";
		}

	} // end all models FOR loop
	gmFile.close();
  
	cout << endl << "Output files: " << pen_file << ".out " << pen_file << ".datasim " <<
		pen_file << ".smod" << endl << endl;
	
	for (uint i=0; i<locus_trees.size(); i++)
		clear_tree(locus_trees[i]);
	locus_trees.clear();
	// clean up
	delete [] GENO_FREQUENCY_LIST;
	GENO_FREQUENCY_LIST = NULL;
	return bestFitness;
}

// ------------------------- F U N C T I O N S --------------------------------

// ----------------------------------------------------------------------------
// This objective function maximizes the difference between table values
// and minimizes the marginal probabilities.
float Objective(GAGenome& g)
{
  GARealGenome& genome = (GARealGenome&) g;
  double marginAvg, heritability, tableVariance, marginVariance, oddsRatio;
  double * marginalProbability = new double[lociCount * NUM_GENOS_PER_LOCUS];

  calc_score(genome, marginalProbability, marginAvg, heritability, tableVariance,
    marginVariance, oddsRatio, CALCODDS, CALCHERIT);
  delete [] marginalProbability;
  
  return calc_fitness(tableVariance, marginAvg, marginVariance, heritability, oddsRatio);
}


// calculate fitness
// Args: Table variance
//       Marginal penetrance average
//       Marginal penetrance variance
//       Heritability
// Ret:  fitness
double calc_fitness(double tableVariance, double marginAvg, double marginVariance,
  double heritability, double oddsRatio){
  
  if(CALCHERIT && isnan(heritability)){
    return -1.0;
  }
  double fitness = 1.0;
      
  // only use table variance when it is set to greater than zero
  if(TABLE_VAR_THRESHOLD > 0 && tableVariance < TABLE_VAR_THRESHOLD) {
    fitness = fitness - (TABLE_VAR_THRESHOLD - tableVariance);
  }
      
  // only use target penetrance when it is set to greater than zero
  if(TARGET_PEN > 0 && fabs(marginAvg - TARGET_PEN) > TARGET_PEN_DIFF){
    fitness = fitness - fabs(marginAvg - TARGET_PEN);  
  }

  // always use marginal penetrance variance in determining fitness
  if(marginVariance > MARGINS_VAR_THRESHOLD) {
    fitness = fitness - (marginVariance - MARGINS_VAR_THRESHOLD)*MARGIN_WEIGHT;
  }
  
  if(CALCHERIT && fabs(heritability - HERITABILITY_THRESHOLD) > .002) {
    fitness = fitness - fabs(heritability - HERITABILITY_THRESHOLD)*HERIT_WEIGHT;
  }  
  
  if(CALCODDS && fabs(oddsRatio - ODDS_TARGET) > .05){
    fitness = fitness - fabs(oddsRatio - ODDS_TARGET) * ODDS_WEIGHT;
  }

  return fitness;

}


// calculate scores
// Args:  genome
//        marginalProbability is the array for storing marginals [out]
//        marginAvg is the average of marginals [out]
//        heritability [out]
//        tableVariance [out]
//        marginVariance [out]
void calc_score(GARealGenome& genome, double * marginalProbability, double & marginAvg,
  double & heritability, double & tableVariance, double & marginVariance, double & oddsRatio,
  bool oddsCalc, bool heritCalc){

  unsigned int genome_length = genome.length();
  unsigned int num_loci = lociCount;
  double * tableValue = new double[genome_length];

  // set table values from genome
  for(unsigned int i=0; i<genome_length; i++)
    tableValue[i] = genome.gene(i);

  if(oddsCalc){
    oddsRatio = calculate_odds_ratio(tableValue, GENO_FREQUENCY_LIST, genome_length);  
  }

  CalculateMarginals(tableValue, marginalProbability, geno_freqs, locus_trees);
  marginVariance = CalculateVariance(marginalProbability, num_loci * NUM_GENOS_PER_LOCUS);
  tableVariance = CalculateVariance(tableValue, genome_length);

  if(heritCalc){
    heritability = CalculateHeritability(locus_trees[0], marginalProbability, 
      num_loci, NUM_GENOS_PER_LOCUS, geno_freqs, tableValue, &marginAvg);
  }
  else if(TARGET_PEN > 0){
    marginAvg = calculate_marginal_mean(geno_freqs, marginalProbability,
      num_loci, NUM_GENOS_PER_LOCUS);
  }
    
  delete [] tableValue;
  
}


// displays statistics to cout as part of interval reports
void show_score_update(GARealGenome& g){
	double marginAvg, heritability, tableVariance, marginVariance, oddsRatio;
	double * marginalProbability = new double[lociCount * NUM_GENOS_PER_LOCUS];
	calc_score(g, marginalProbability, marginAvg, heritability, tableVariance,
		marginVariance, oddsRatio, true, true);

	cout << setw(14) <<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << marginVariance;
	cout << setw(14) <<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << tableVariance;
	cout << setw(14) <<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << heritability;
	cout << setw(12) <<setiosflags(ios::fixed |ios::showpoint)<<setprecision(6) << marginAvg;
	cout << setw(12) <<setiosflags(ios::scientific|ios::showpoint)<<setprecision(6) << oddsRatio;
}


// calculates marginal means and heritability on the genome passed in
// used for showing scores on best individuals
void summarize_scores(GARealGenome& genome, ofstream & best_model_file, SimPreferences & pref,
	vector< vector< tree<node *> > >  & subModelTrees){
	double marginAvg, heritability, tableVariance, marginVariance, oddsRatio;
	double * marginalProbability = new double[lociCount * NUM_GENOS_PER_LOCUS];
	
	calc_score(genome, marginalProbability, marginAvg, heritability, tableVariance,
		marginVariance, oddsRatio, true, true);
    
	if(pref.get_verbose()){
		best_model_file << "Marginals = " ;
		for(unsigned int i=0; i<lociCount * NUM_GENOS_PER_LOCUS; i++){
			best_model_file << marginalProbability[i] << " " ;
		}
		best_model_file << endl << endl;
	} 
  
	best_model_file << setw(15) << setfill(' ') << "Model" << setw(15) << "MarginalVar";
	if(pref.get_verbose()){
    	best_model_file << setw(12) << "TableVar" << setw(14) << "Heritability"
      		<< setw(12) << "MarginAvg" << setw(12) << "Odds Ratio" << setw(10) << "Fitness";
	}
	best_model_file << endl;
	best_model_file << "---------------------------------------------------------------------------------------" << endl;    

	// show loci in model  
	char first_locus = 'A';
	string lociModel;
	for(unsigned int z=0; z<lociCount; z++){
		lociModel += char(first_locus + z);
		lociModel += " ";
	}
	best_model_file << setw(15) << lociModel << setw(15) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << marginVariance;

	if(pref.get_verbose()){
		double fitness = calc_fitness(tableVariance, marginAvg, marginVariance, heritability, oddsRatio);
    	best_model_file << setw(12) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(8) << tableVariance;
		best_model_file << setw(14) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(8) << heritability;
    	best_model_file << setw(12) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << marginAvg;
		best_model_file << setw(12) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << oddsRatio;
		best_model_file << setw(10) << setiosflags(ios::fixed|ios::showpoint)<<setprecision(6) << fitness;    
	}
	best_model_file << endl;
	
	//check for submodels
	if(pref.get_submodels_on() && lociCount > 2){
		unsigned int genome_length = genome.length();
		double * tableValue = new double[genome_length];
		for(unsigned int i=0; i<genome_length; i++)
			tableValue[i] = genome.gene(i);
		SubModelCheck(tableValue, genome_length, best_model_file, pref, subModelTrees);  
		delete [] tableValue;
	}
	
	delete [] marginalProbability;
	best_model_file << endl << endl;;  
}


void SubModelCheck(double * tableValue, unsigned int genome_length,  ofstream & gmFile, 
  SimPreferences & pref, vector< vector< tree<node *> > >  & subModelTrees){

  vector<int> lociList;
  for(unsigned int i=0; i<lociCount; i++){
    lociList.push_back(i);
  }

  map<string, resultInfo> resultMap;

  checkSubModels(tableValue, lociCount, lociList, 
    resultMap, geno_freqs, subModelTrees);

  // output submodel loci
  char first_locus = 'A';
  
  map<string, resultInfo>::iterator mapIter;

  for(mapIter=resultMap.begin(); mapIter != resultMap.end(); mapIter++){
    int numLoci = (mapIter->second).lociList.size();
    string lociCombination;
    for(int i=0; i<numLoci; i++){
      lociCombination += char(first_locus + (mapIter->second).lociList[i]);
      lociCombination += " ";
    }
    gmFile << setw(15) << lociCombination << setw(15) << (mapIter->second).marginVariance;
    if(pref.get_verbose()){
      gmFile << setw(12) << (mapIter->second).tableVariance << setw(14) 
        << (mapIter->second).herit;
      double fitness = calc_fitness((mapIter->second).tableVariance, (mapIter->second).marginAvg, 
        (mapIter->second).marginVariance, (mapIter->second).herit, (mapIter->second).oddsRatio);      
      gmFile << setw(12) << (mapIter->second).marginAvg << setw(12) << (mapIter->second).oddsRatio
      << setw(10) << fitness;        
    }
    gmFile << endl;
    // To prevent memory leak need to delete dynamically allocated marginal penetrances array
    delete [] (mapIter->second).margPenetrances;
  }
  gmFile << endl;
}


void AddModelToBest(ofstream & best_model_file, GARealGenome& genome, int num_loci,
   const vector<string> & labels){
   int i;

   best_model_file.flags(ios::left|ios::showpoint);
   int originalPrecision = best_model_file.precision();
   for(i=0; i < genome.length(); i++) {
      if(genome.gene(i) < .01)
        best_model_file.precision(1);
      else if(genome.gene(i) < .1)
        best_model_file.precision(2);
      else
        best_model_file.precision(3);
        best_model_file <<  labels[i] << " " << setw(5) << setfill('0') << genome.gene(i) << endl;
   }
   best_model_file << endl;
   best_model_file.precision(originalPrecision);
   best_model_file.flags(ios::left|ios::showpoint);
}


// function: CreateDatasimFile  -- creates sample penetrance file with penetrances
//           of model for use in data simulator program
// args:  name of datasim file
//        genome
//        number of loci
//        allele frequencies
// ret:   name of datasim file
string CreateDatasimFile(string base_file, GARealGenome& genome, int num_loci, 
  const vector< vector<double> > & all_freqs, SimPreferences & pref, int modNum){
	stringstream datasim_file;
	datasim_file<<base_file<<"."<<modNum<<".datasim";
  ofstream datasimFile(datasim_file.str().c_str(), ios::out);

  datasimFile << "# this file is needed for library data simulation config file" << std::endl << std::endl;
  datasimFile << "# seed for random number generator" << std::endl;
  datasimFile << "RAND " << pref.get_seed() << std::endl << std::endl;
  datasimFile << "# Each model needs its own file showing the penetrance table" << std::endl;
  datasimFile << "# use keyword and then number of model files after" << std::endl;
  datasimFile << "# followed by fraction of individuals determined using this" << std::endl;
  datasimFile << "# model (total should be 1.0)." << std::endl;
  datasimFile << "# If you want all status to be randomly assigned," << std::endl;
  datasimFile << "# remove or comment out all the MODELFILES and" << std::endl;
  datasimFile << "# set MODELFILES to 0." << std::endl;
  datasimFile << "MODELFILES 1" << std::endl;
  datasimFile <<  base_file << "." << modNum << ".smod 1.0" << std::endl << std::endl;
  datasimFile << "# number of datasets created using this file" << std::endl;
  if(pref.get_simsets() > 0)
    datasimFile << "SIMSETS " << pref.get_simsets() <<  std::endl << std::endl;
  else
    datasimFile << "SIMSETS 1" <<  std::endl << std::endl;
  datasimFile << "# percentage of error in genotypes" << std::endl;
  datasimFile << "GENOTYPEERROR " << pref.get_genoerror() << std::endl << std::endl;
  datasimFile << "# fraction of affected that are due to environmental factors" << std::endl;
  datasimFile << "PHENOCOPY " << pref.get_phenocopy() << std::endl << std::endl;
  
  if(pref.get_simtype() == ProbabilitySimType){
    datasimFile << "# simulation type" << std::endl;
    datasimFile << "SIMTYPE prob" << std::endl << std::endl; 
    datasimFile << "# number of affected to simulate" << std::endl;
    datasimFile << "AFFECTED " << pref.get_num_affected() << std::endl << std::endl;
    datasimFile << "# unaffected to simulate" << std::endl;
    datasimFile << "UNAFFECTED " << pref.get_num_unaffected() << std::endl << std::endl;
    datasimFile << "# total number of loci to simulate" << std::endl;
    datasimFile << "SIMLOCI " << pref.get_simloci() << std::endl << std::endl;
  }
  else if(pref.get_simtype() == PopulationSimType ){
    datasimFile << "# simulation type" << std::endl;
    datasimFile << "SIMTYPE pop" << std::endl << std::endl;
    datasimFile << "# population size" << std::endl;
    datasimFile << "POPSIZE " << pref.get_simpopsize() << std::endl << std::endl;;
    datasimFile << "# number of generations to simulate" << std::endl;
    datasimFile << "NUMGENS " << pref.get_numgens() << std::endl << std::endl;;
    datasimFile << "# number of generations to simulate" << std::endl;
    datasimFile << "GENES " << pref.get_numgenes() << std::endl << std::endl;
    datasimFile << "# minimum number of SNPs per gene" << std::endl;
    datasimFile << "# if want every gene to have same number of SNPs make MINSNP and MAXSNP equal" << std::endl;
    datasimFile << "MINSNP " << pref.get_minsnp() << std::endl << std::endl;
    datasimFile << "# maximum number of SNPs per gene" << std::endl;
    datasimFile << "MAXSNP " << pref.get_maxsnp() << std::endl << std::endl;  
    datasimFile << "# minimum recombination rate between SNPS in a gene" << std::endl;
    datasimFile << "MINRECOMB " << pref.get_minrecomb() << std::endl << std::endl;
    datasimFile << "# maximum recombination rate between SNPs in a gene" << std::endl;
    datasimFile << "MAXRECOMB " << pref.get_maxrecomb() << std::endl << std::endl;
  }
  
  datasimFile << "# For setting random allele frequencies" << std::endl;
  datasimFile << "# list minor allele frequency minimum and maximum" << std::endl;
  datasimFile << "# when this isn't set then use DEFAULTALLELE to make" << std::endl;
  datasimFile << "# every allele that isn't listed under ALLELEFREQS the same" << std::endl;
  std::vector<double> allelelimits = pref.get_allelelimits();
  if(allelelimits.size() > 0)
    datasimFile << "ALLELELIMITS " << allelelimits[0] << " " << allelelimits[1]<< std::endl;

  datasimFile << std::endl << "# Default allele Frequencies (use if not using random allele freqs)" << std::endl;
  datasimFile << "# this only applies if ALLELELIMITS is commented out or missing " << std::endl;
  datasimFile << "# in the file" << std::endl;
  datasimFile << "DEFAULTALLELE " << all_freqs[0][0] << " " << all_freqs[0][1] << std::endl << std::endl;
  datasimFile << "# for specified allele frequencies " << std::endl;
  datasimFile << "# first is locus (start with 0) followed by major and minor" << std::endl;
  datasimFile << "# allele frequencies" << std::endl;
  datasimFile << "# should specify frequencies for the disease allele" << std::endl;
  datasimFile << "# loci in model files" << std::endl;
  datasimFile << "ALLELEFREQS" << std::endl;
  int startlocus = 5;
  for(int i=0; i<num_loci; i++)
    datasimFile << startlocus + i*5 << " " << all_freqs[i][0] << " " << all_freqs[i][1] << std::endl;  
  datasimFile.close();
  return datasim_file.str();
}

//function:  CreateModelFile -- currently assumes all loci have same frequencies
//           For use with new data simulation library
//args:  numloci
//    :  vector< vector <double> > geno_freqs
string CreateModelFile(string base_file, int num_loci, GARealGenome& genome,
  const vector<string> & labels, int modNum, ModelLociArray &loci){
	stringstream freq_file;
	freq_file<<base_file<<"."<<modNum<<".smod";
	ofstream fFile(freq_file.str().c_str(), ios::out);
	
	fFile << "# specify which loci are disease loci" << std::endl;
	fFile << "DISEASELOCI ";
	// add loci spaced 5 apart
	
	char al1 = 'A';
	char al2 = 'a';
 
	fFile<<"# This penetrance table was generated by simPEN.\n";
	fFile<<"FREQ_THRESHOLD 0.00001\n\n";
  	for(int i=0; i<num_loci; i++){
		fFile<<"# Locus "<<i<<", "<<loci[i].label<<"\n";
		fFile<<"FREQ "<<al1++<<" "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<loci[i].alFreq1<<"\n";
		fFile<<"FREQ "<<al2++<<" "<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)<<loci[i].alFreq2<<"\n";
   	}
	fFile << std::endl << std::endl;
	
	fFile << "PENTABLE" << std::endl; 
	AddModelToBest(fFile, genome, num_loci, labels);
	
	fFile.close();
	return freq_file.str();
}


// function: CreateLabels
//           Labels are in order of AABB, AABb, AAbb with the outermost 
//           loci always cycling the fastest
// Arg:  Number of loci
//       vector reference to add labels
// Ret:  none
void CreateLabels( int num_loci, vector<string> & labels){
   char upper_start = 'A';
   char lower_start = 'a';

   vector< vector<string> > genotypes;
   int i;
   unsigned int * lower_index = new unsigned int[num_loci];
   unsigned int * upper_index = new unsigned int[num_loci];

   for(i=0; i<num_loci; i++){
      genotypes.push_back(vector<string>(3));
      string upper(1,upper_start);
      string lower(1,lower_start);
      genotypes[i][0] = upper + upper;
      genotypes[i][1] = upper + lower;
      genotypes[i][2] = lower + lower;
      upper_start++;
      lower_start++;
      lower_index[i] = 0;
      upper_index[i] = 3;
   }

   header_nested_loops(lower_index, upper_index, num_loci, labels, genotypes);

   delete [] lower_index;
   delete [] upper_index;
}


void header_nested_loops( unsigned int * lower, unsigned int *upper, int max_depth,
     vector<string> & labels, vector< vector<string> > & genotypes){

   unsigned int * indexes = new unsigned int[max_depth];
   int cur_depth = 0;

   indexes[cur_depth] = lower[cur_depth];

   while(1){
      if( indexes[cur_depth] < upper[cur_depth] ) {
         if( cur_depth == max_depth - 1 ) {
            string new_label;
           // create label
           for(int i =0; i<=cur_depth; i++){
                new_label += genotypes[i][indexes[i]];
            }
            labels.push_back(new_label);
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

// ----------------------------------------------------------------------------

// Outputs information on program to user
void header(string progname, string versiondate ){
        cout << endl;
        cout << "---------------------------------------------------" << endl;
        cout << "       " << "Name:    " << progname <<  endl;
        cout << "       " << "Date:    " << versiondate <<  endl;
        cout << "       " << "Usage:   " << progname << " <config file> [output name]"<< endl << endl;
        cout << "       " << "Example: " << progname << " sample.simpen 2locus" << endl;
        cout << "---------------------------------------------------" << endl << endl;
}

}


//

