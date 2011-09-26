// SimPreferences.h

// Reads preferences file and sets the preferences

#ifndef __SIMPREFERENCES_H__
#define __SIMPREFERENCES_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
//#include "simpen/diseaselocus.h"
#include "simulation/diseaselocus.h"

namespace SimPen {


using namespace std;


const string LOCI = "LOCI", SIMPENSEED = "SEED", FREQ = "FREQ", TABLEVAR = "TABLEVAR",
   MARGVAR = "MARGVAR", HERIT = "HERIT", MODELS = "MODELS", GEN = "GEN",
   POPSIZE = "POPSIZE", DEMES = "DEMES", MUTATE = "MUTATE", CROSS = "CROSS",
   VERBOSE = "VERBOSE", UPDATE = "UPDATE", USETABLEVAR = "USETABLEVAR", PENTARGET = "PENTARGET",
   SUBMODELS = "SUBMODELS", MARGWEIGHT = "MARGWEIGHT", HERITWEIGHT = "HERITWEIGHT",
   SIMPENNUMSIMSETS = "SIMSETS", PHENOCOPYRATE = "PHENOCOPY", SIMGENOERROR = "GENOERROR",
   SIMDATAFILEBASE = "SIMBASE", SIMAFFECTED = "AFFECTED", SIMUNAFFECTED = "UNAFFECTED", SIMDATALOCI = "SIMLOCI",
   LOCIPOS = "LOCIPOS" , SIMDATASIMFILE = "MODELFILES", SIMALLELELIMITS = "ALLELELIMITS",
   GENOSIMTYPE = "SIMTYPE", SIMPOPSIZE = "SIMPOPSIZE", SIMNUMGENS = "NUMGENS", SIMGENES = "GENES",
   SIMMINSNP = "MINSNP", SIMMAXSNP = "MAXSNP", SIMMINRECOMB = "MINRECOMB", 
   SIMMAXRECOMB = "MAXRECOMB", ODDSRATIO = "ODDSRATIO", ODDSWEIGHT = "ODDSWEIGHT", ADDLOCUS="ADD_LOCUS";

const string DEFAULT_VERBOSE = "OFF", DEFAULT_USETABLEVAR= "OFF", DEFAULT_SUBMODELS = "OFF";
const string DEFAULT_DATAFILEBASE = "simpen";
const string DEFAULT_SIMDATASIMFILE = "";
const string DEFAULT_SIMTYPE = "prob";

enum genoSimTypes{ 
  ProbabilitySimType,
  PopulationSimType
};
    

class SimPreferences{

  public:
    SimPreferences();
    SimPreferences(string pref_file);

    // Use to read in configuration file
    void read_config(string configname); // throws PrefExcept

    // Returns seed number
    unsigned int get_seed() const;

    //Returns number of loci
    //int get_loci() const;

    // Returns 2-D vector of allele frequencies
    //vector< vector<double> > get_freqs();

    int get_num_models() const;
    unsigned int get_gen()const;
    unsigned int get_demes()const;
    unsigned int get_pop_size()const;

    double get_target_penetrance()const;
    double get_marg_var()const;
    double get_herit()const;
    double get_mutate()const;
    double get_cross()const;
    double get_table_var()const;
    bool get_verbose()const{return verbose_on;}
    bool get_tablevar_on()const{return tablevar_on;}
    int get_update_interval()const;
    bool get_submodels_on()const{return submodels_on;}
    int get_marg_weight()const{return marginal_weight;}
    int get_herit_weight()const{return herit_weight;}
    int get_num_affected()const{return num_affected;}
    int get_num_unaffected()const{return num_unaffected;}
    double get_phenocopy() const{return phenocopy;}
    double get_genoerror() const{return genoerror;}

    int get_simpopsize()const{return simpopulation;}
    int get_numgens() const{return numgenerations;}
    int get_numgenes() const{return numgenes;}
    int get_minsnp() const{return minnumsnp;}
    int get_maxsnp() const{return maxnumsnp;}
    double get_minrecomb() const{return minrecomb;}
    double get_maxrecomb() const{return maxrecomb;}

    bool datasimfile_set() const;
    string get_datasimfile() const{return simDataSimFile;}
    int get_simloci() const{return simloci;}
    int get_simsets() const{return num_simsets;}
    double get_odds_ratio() const {return odds_ratio;}
    double get_odds_weight() const {return odds_weight;}
    genoSimTypes get_simtype() const{return simulationType;}
    std::vector<double> get_allelelimits() const{return simAlleleLimits;}
    std::vector<int> get_simloc_pos() const {return simlocpos;}

	//void set_loci(int lociCount) ;
//	void AppendFrequencies(float freq1, float freq2);
		
	void SetSeed(int seed) { rand_seed=seed; }
	
	//void ResetAlleleFreq();

	/**
	 * @brief Add a new locus to the disease model
	 */
	void AddDiseaseLocus(const char *, int chr, int idx, float fr1, float fr2);

	void AddDiseaseLocus(Simulation::StatusModel::DiseaseLocus &l);

	/**
	 * @brief Return a locus from the current model
	 */
	Simulation::StatusModel::DiseaseLocus &GetDiseaseLocus(int idx);

	/**
	 * @brief Returns the number of loci in the current model	
 	 */
	int GetLociCount();

	/**
	 * @brief Removes all elements of the model
	 */
	void ResetModel();
  private:
    // assigns values to parameters based on reading of configuration file
    void set_params();

    // skips remainder of a line after parameters are read
    void skip_rest(ifstream & filestream);

    // returns the value as a string of the parameter in the log file
    string get_param(string key, map<string, string> & params);

    // sets the variables to be equal to values in config file
    void set_params(map<string, string> & pref_map);

    // checks that parameters are set and within range
    void check_params();

    // returns true when string passed in equals ON
    bool check_on_off(string status);

    // converts long to string
    string itos(long number);
    string itos(float number);
    string itos(int number);
    string itos(unsigned int number);
    string itos(double number);

    unsigned int rand_seed, init_seed,num_gens, num_demes, pop_size;
    double herit, marg_var, table_var, mutation_rate, crossover_rate, target_pen,
      phenocopy, genoerror, minrecomb, maxrecomb, odds_ratio, odds_weight;
	//int num_loci;
    int num_models, update_interval, marginal_weight, herit_weight,
      num_simsets, num_affected, num_unaffected, simloci, simpopulation,
      numgenerations, numgenes, minnumsnp, maxnumsnp;
    string pref_file, simDataSimFile;
    bool verbose_on, tablevar_on, submodels_on;

	Simulation::StatusModel::ModelLociArray diseaseLoci;

    //std::vector< std::vector<double> > allele_freqs;
    std::vector<double> simAlleleLimits;
    std::vector<int> simlocpos;
    genoSimTypes simulationType;
    map<string, genoSimTypes> genoSimTypesMap;


   	static int const NUMLOCI = 2;
    static unsigned int const SEED_MAX = INT_MAX;
    static int const MIN_LOCI = 2;
    // set up allele frequencies
#define DEFAULT_FREQ_HIGH 0.8
#define DEFAULT_FREQ_LOW 0.2
    static int const DEFAULT_POP = 100;
    static int const DEFAULT_DEMES = 10;
    static int const DEFAULT_GENS = 5000;
    int static const DEFAULT_MODELS = 1;
    int static const DEFAULT_UPDATE = 100;
#define DEFAULT_MARG 0.00001
#define DEFAULT_HERIT 0.0
#define DEFAULT_TABLE 0.1
//#define DEFAULT_MAXSNP 5
//#define DEFAULT_MINSNP 5
#define DEFAULT_MUTATION 0.001
#define DEFAULT_CROSSOVER 0.9
#define DEFAULT_GENOERROR 0.0

#define DEFAULT_MINRECOMB 0.005
#define DEFAULT_MAXRECOMB 0.005
#define DEFAULT_PHENORATE 0.0
#define DEFAULT_TARGETPEN 0.0
#define DEFAULT_ODDS_RATIO 0
#define DEFAULT_ODDS_WEIGHT 5


    int static const DEFAULT_MARGWEIGHT = 5;
    int static const DEFAULT_HERITWEIGHT = 5;
    
    int static const DEFAULT_SIMSETS = 0;
    int static const DEFAULT_AFFECTED = 200;
    int static const DEFAULT_UNAFFECTED = 200;
    int static const DEFAULT_SIMLOCI = 20;
    
    int static const DEFAULT_MINSNP = 5;
    int static const DEFAULT_MAXSNP = 5;
    int static const DEFAULT_SIMPOP = 200;
    int static const DEFAULT_NUMGENES = 10;
    int static const DEFAULT_NUMGENERATIONS = 10;
    
    
};

}

#endif
