// SimPreferences.cpp

#include "prefexcept.h"

#include "simpreferences.h"
#include <math.h>

#ifdef WIN32
#include <windows.h>
#include <sys/timeb.h>
#include <time.h>
#endif

namespace SimPen {

using namespace Simulation::StatusModel;

SimPreferences::SimPreferences(){
  init_seed = (unsigned int)(time(NULL));
}

SimPreferences::SimPreferences(string pref_file){
  init_seed = (unsigned int)(time(NULL));
  read_config(pref_file);
}

 // Returns seed number
unsigned int SimPreferences::get_seed() const{
  return rand_seed;
}

/*Returns number of clusters
int SimPreferences::get_loci() const{
  return num_loci;
}

void SimPreferences::set_loci(int lociCount) {	
	num_loci = lociCount;
}

 // Return allele frequencies
 vector< vector<double> > SimPreferences::get_freqs(){
   return allele_freqs;
 }

void SimPreferences::ResetAlleleFreq() {
	allele_freqs.clear();
}

void SimPreferences::AppendFrequencies(float freq1, float freq2) {
	vector<double> af;
	af.push_back(freq1);
	af.push_back(freq2);
	allele_freqs.push_back(af);
}*/

void SimPreferences::ResetModel() {
	diseaseLoci.clear();
}

void SimPreferences::AddDiseaseLocus(DiseaseLocus &l) {
	diseaseLoci.push_back(l);
}

void SimPreferences::AddDiseaseLocus(const char *lbl, int chr, int idx, float fr1, float fr2) {
	diseaseLoci.push_back(DiseaseLocus(lbl, chr, idx, fr1, fr2));
}

DiseaseLocus &SimPreferences::GetDiseaseLocus(int idx) {
	//assert((uint)idx < diseaseLoci.size());
	return diseaseLoci[idx];
}

int SimPreferences::GetLociCount() {
	return diseaseLoci.size();
}

// returns a map containing the parameters read in from the config
// file
void SimPreferences::read_config(string configname = ""){
  map<string, string> pref_map;
  
  //set SimTypesMap
  genoSimTypesMap["pop"] = PopulationSimType;
  genoSimTypesMap["prob"] = ProbabilitySimType;

  //set defaults
  pref_map[SIMPENSEED] = itos(init_seed);
  pref_map[LOCI] = itos(NUMLOCI);
  pref_map[TABLEVAR] = itos(DEFAULT_TABLE);
  pref_map[MARGVAR] = itos(DEFAULT_MARG);
  pref_map[HERIT] = itos(DEFAULT_HERIT);
  pref_map[MODELS] = itos(DEFAULT_MODELS);
  pref_map[GEN] = itos(DEFAULT_GENS);
  pref_map[POPSIZE] = itos(DEFAULT_POP);
  pref_map[DEMES] = itos(DEFAULT_DEMES);
  pref_map[MUTATE] = itos(DEFAULT_MUTATION);
  pref_map[CROSS] = itos(DEFAULT_CROSSOVER);
  pref_map[UPDATE] = itos(DEFAULT_UPDATE);
  pref_map[PENTARGET] = itos(DEFAULT_TARGETPEN);
  
  pref_map[VERBOSE] = DEFAULT_VERBOSE;
  pref_map[USETABLEVAR] = DEFAULT_USETABLEVAR;
  pref_map[SUBMODELS] = DEFAULT_SUBMODELS;
  
  pref_map[MARGWEIGHT] = itos(DEFAULT_MARGWEIGHT);
  pref_map[HERITWEIGHT] = itos(DEFAULT_HERITWEIGHT);

  pref_map[SIMPENNUMSIMSETS] = itos(DEFAULT_SIMSETS);
  pref_map[PHENOCOPYRATE] = itos(DEFAULT_PHENORATE);
  pref_map[SIMGENOERROR] = itos(DEFAULT_GENOERROR);
  pref_map[SIMUNAFFECTED] = itos(DEFAULT_UNAFFECTED);
  pref_map[SIMDATALOCI] = itos(DEFAULT_SIMLOCI);
  pref_map[SIMDATASIMFILE] = DEFAULT_SIMDATASIMFILE;
  
  pref_map[SIMMINSNP] = itos(DEFAULT_MINSNP);
  pref_map[SIMMAXSNP] = itos(DEFAULT_MAXSNP);
  pref_map[SIMMINRECOMB] = itos(DEFAULT_MINRECOMB);
  pref_map[SIMMAXRECOMB] = itos(DEFAULT_MAXRECOMB);
  pref_map[SIMGENES] = itos(DEFAULT_NUMGENES);
  pref_map[SIMPOPSIZE] = itos(DEFAULT_SIMPOP);
  pref_map[SIMNUMGENS] = itos(DEFAULT_NUMGENERATIONS);
  
  pref_map[ODDSRATIO] = itos(DEFAULT_ODDS_RATIO);
  pref_map[ODDSWEIGHT] = itos(DEFAULT_ODDS_WEIGHT);

  ifstream config(configname.c_str(), ios::in);
  if(!config.is_open()){
    throw PrefExcept(configname + ":  cannot open!");
  }

  string key, value;

  // read in each parameter to use as key in map and its parameter
  while(!config.eof()){
    config >> key;
    
    // make sure keyword is uppercase
    for(unsigned int strpos=0; strpos<key.length(); strpos++){
      key[strpos] = toupper(key[strpos]);
    }
    
	if (key.compare(ADDLOCUS) == 0) {
		int chr, idx;
		float af1, af2;
		
		config>>chr>>idx>>af1>>af2;
		stringstream ss;
		ss<<chr<<":"<<idx;
		AddDiseaseLocus(ss.str().c_str(), chr, idx, af1, af2);
	}
   /*if(key.compare(FREQ) == 0){
      allele_freqs.push_back(vector<double>(2,0));
      config >> allele_freqs[loci_freq][0] >> allele_freqs[loci_freq][1];
      loci_freq++;
    }*/
    else if(key.compare(LOCIPOS) == 0){
      char line[512];
      config.getline(line, 512);
      string split(line);
      
      int lastDig = split.find_last_of("0123456789");
      split = split.substr(0, lastDig+1);
      stringstream ss(split);
      while(!ss.eof()){
        int position;
        ss >> position;
        simlocpos.push_back(position);
      }
      
      
    }
    else if(key.compare(SIMALLELELIMITS) ==0){
      char line[512];
      config.getline(line, 512);
      string split(line);
      int lastDig = split.find_last_of("0123456789");
      split = split.substr(0, lastDig+1);
      stringstream ss(split);
      while(!ss.eof()){
        double allLimit;
        ss >> allLimit;
        simAlleleLimits.push_back(allLimit);
      }
    }
    else{
      config >> value;
      pref_map[key] = value; // insert pair into map 
    }
    skip_rest(config); // skip any comments after values
    key="";
  }

  config.close();
  if(pref_map.find(GENOSIMTYPE) == pref_map.end()){
    if(pref_map.find(SIMAFFECTED) == pref_map.end()){
      pref_map[GENOSIMTYPE] = "pop";
    }
    else{
      pref_map[GENOSIMTYPE] = "prob";
    }
  }

  if(pref_map.find(SIMAFFECTED) == pref_map.end()){
    pref_map[SIMAFFECTED] = itos(DEFAULT_AFFECTED);
  }

  // now set the parameters using the map
  set_params(pref_map);
}


// sets variables within object to match parameters
// throws LapExcept if no parameter found
void SimPreferences::set_params(map<string, string> & pref_map){
        // construct a string to use in assigning values
        string par_holder;

//   par_holder = get_param(LOCI, pref_map) + ' ';
   par_holder = get_param(SIMPENSEED, pref_map) + ' ';
   par_holder += get_param(TABLEVAR, pref_map) + ' ';
   par_holder += get_param(MARGVAR, pref_map) + ' ';
   par_holder += get_param(HERIT, pref_map) + ' ';
   //par_holder += get_param(MODELS, pref_map) + ' ';
   par_holder += get_param(GEN, pref_map) + ' ';
   par_holder += get_param(POPSIZE, pref_map) + ' ';
   par_holder += get_param(DEMES, pref_map) + ' ';
   par_holder += get_param(MUTATE, pref_map) + ' ';
   par_holder += get_param(CROSS, pref_map) + ' ';
   par_holder += get_param(UPDATE, pref_map) + ' ';
   par_holder += get_param(PENTARGET, pref_map) + ' ';
   par_holder += get_param(MARGWEIGHT, pref_map) + ' ';
   par_holder += get_param(HERITWEIGHT, pref_map) + ' ';
   par_holder += get_param(SIMPENNUMSIMSETS, pref_map) + ' ';
   par_holder += get_param(PHENOCOPYRATE, pref_map) + ' ';
   par_holder += get_param(SIMGENOERROR, pref_map) + ' ';
   par_holder += get_param(SIMAFFECTED, pref_map) + ' ';
   par_holder += get_param(SIMUNAFFECTED, pref_map) + ' ';
   par_holder += get_param(SIMDATALOCI, pref_map) + ' ';

   par_holder += get_param(SIMPOPSIZE, pref_map) + ' ';
   par_holder += get_param(SIMNUMGENS, pref_map) + ' ';
   par_holder += get_param(SIMGENES, pref_map) + ' ';
   par_holder += get_param(SIMMINSNP, pref_map) + ' ';
   par_holder += get_param(SIMMAXSNP, pref_map) + ' ';
   par_holder += get_param(SIMMINRECOMB, pref_map) + ' ';
   par_holder += get_param(SIMMAXRECOMB, pref_map) + ' ';
   par_holder += get_param(ODDSRATIO, pref_map) + ' ';
   par_holder += get_param(ODDSWEIGHT, pref_map) + ' ';
   par_holder += get_param(SIMDATASIMFILE, pref_map) + ' ';

   
  stringstream ss(par_holder);

	//Saving more than one model / run really doesn't mean much any more
	//ss >> num_loci;
	num_models=1;
  ss >> rand_seed >> table_var >> marg_var >> herit;// >> num_models
  ss >> num_gens >> pop_size >> num_demes >> mutation_rate >> crossover_rate 
   >> update_interval >> target_pen >> marginal_weight >> herit_weight 
   >> num_simsets >> phenocopy >> genoerror >> num_affected >> num_unaffected
   >> simloci >> simpopulation >> numgenerations >> numgenes 
   >> minnumsnp >> maxnumsnp >> minrecomb >> maxrecomb 
   >> odds_ratio >> odds_weight >> simDataSimFile;

  verbose_on = check_on_off(get_param(VERBOSE, pref_map));
  tablevar_on = check_on_off(get_param(USETABLEVAR, pref_map));
  submodels_on = check_on_off(get_param(SUBMODELS, pref_map));

  /* when don't have enough frequencies passed in repeat the first
  // one if there is one, else use the default frequencies
  double low_freq = 0.0;
  double high_freq = 0.0;
  if(allele_freqs.size() > 0){
   high_freq = allele_freqs[0][0];
   low_freq = allele_freqs[0][1];
  }
  else{
   low_freq = DEFAULT_FREQ_LOW;
   high_freq = DEFAULT_FREQ_HIGH;
  }
  int count = allele_freqs.size();
	
  if(num_loci < NUMLOCI)
    throw PrefExcept(LOCI + " must be set to at least 2");

  while(count != num_loci){
    allele_freqs.push_back(vector<double>(2,0));
    allele_freqs[count][0] = high_freq;
    allele_freqs[count][1] = low_freq;
    count++;
  }*/

  // when loci pos not set use default
  if(simlocpos.size() ==0){
    int startloc = 4;
	int lociCount = GetLociCount();
    for(int currLoc=0; currLoc < lociCount; currLoc++){
      simlocpos.push_back(startloc + currLoc*5);
    }
  }
  
  // set simtype
  // check that GENOSIMTYPE is acceptable 
  map<string, genoSimTypes>::iterator simtypeIter;

  simtypeIter = genoSimTypesMap.find(get_param(GENOSIMTYPE, pref_map));
  if(simtypeIter != genoSimTypesMap.end())
    simulationType = simtypeIter->second;
  else
    throw PrefExcept(GENOSIMTYPE + " must be either pop or prob");
  
  if(!datasimfile_set()){
    check_params();
  }
}


void  SimPreferences::check_params(){
/* EST Doesn't work for genomeSIMLA
  if(GetLociCount() < 2)
    throw PrefExcept("A disease model requires at least 2 loci to be epistatic");
*/
  if(rand_seed > SEED_MAX || rand_seed < 0)
    throw PrefExcept(SIMPENSEED + " must be greater than 0 and less than " + itos(SEED_MAX));
  if(table_var <= 0 || table_var > 1)
   throw PrefExcept(TABLEVAR + " must be greater than 0 and less than 1");
  if(marg_var <= 0 || marg_var > 1)
   throw PrefExcept(MARGVAR + " must be greater than 0 and less than 1");
  if(herit < 0 || herit > 1)
   throw PrefExcept(HERIT + " must be greater than or equal to 0 and less than 1");
  if(crossover_rate < 0 || crossover_rate > 1)
   throw PrefExcept(CROSS + " must be greater than 0 and less than 1");
  if(mutation_rate < 0 || mutation_rate > 1)
   throw PrefExcept(MUTATE + " must be greater than or equal to 0 and less than 1");
  if(num_models <= 0)
   throw PrefExcept(MODELS + " must be greater than 0");
//  if(num_loci <= 0)
//   throw PrefExcept(LOCI + " must be greater than 0");
  if(odds_ratio < 0)
   throw PrefExcept (ODDSRATIO + " must be greater than 0");
  if(herit_weight < 1 || herit_weight > 1000000)
   throw PrefExcept(HERITWEIGHT + " must be between 1 and 1000000");
  if(marginal_weight < 1 || marginal_weight > 1000000)
   throw PrefExcept(MARGWEIGHT + " must be between 1 and 1000000");
  if(odds_weight < 0 || odds_weight > 1000000)
   throw PrefExcept(ODDSWEIGHT + " must be between 0 and 1000000");
  if(num_demes <= 0)
   throw PrefExcept(DEMES + " must be greater than 0");
  if(num_gens <= 0)
   throw PrefExcept(GEN + " must be greater than 0");
  if(pop_size <= 0)
   throw PrefExcept(POPSIZE + " must be greater than 0");
  if(verbose_on !=0 && verbose_on != 1)
   throw PrefExcept(VERBOSE + " must be either 0 (off) or 1 (on)");
  if(update_interval <= 0)
   throw PrefExcept(UPDATE + " must be greater than 0");
   // check frequencies
	int lociCount = GetLociCount();
	for (int locus=0; locus<lociCount; locus++) {
		if (!diseaseLoci[locus].Verify()) {
			diseaseLoci[locus].Report();
			throw PrefExcept("Allele Frequency for Loci must be equal to one");
		}
	}
   /*
	for(int locus=0; locus < num_loci; locus++)
      if(fabs(1.0 - (allele_freqs[locus][0] + allele_freqs[locus][1])) > 0.001)
        throw PrefExcept("All " + FREQ + " lines must total 1.0 for the two frequencies");
	*/

/* This is all simulation stuff which isn't particularly important

  if(phenocopy <0 || phenocopy > 1)
    throw PrefExcept(PHENOCOPYRATE + " must be equal to or greater than 0 and less than or equal to 1.0");
  if(genoerror < 0 || genoerror > 1)
    throw PrefExcept(SIMGENOERROR + " must be equal to or greater than 0 and less than or equal to 1.0");
  if(num_affected < 0)
    throw PrefExcept(SIMAFFECTED + " must be equal to or greater than 0");
  if(num_unaffected < 0)
    throw PrefExcept(SIMUNAFFECTED + " must be equal to or greater than 0");
  if(simloci < num_loci){
    throw PrefExcept(SIMDATALOCI + " must be equal to or greater than number of loci in model "+
      itos(num_loci));
  }
  
  if(num_simsets < 0)
    throw PrefExcept(SIMPENNUMSIMSETS + " must be greater than 0");
*/  
  // check that all entries for position of loci are within possible range
  for(unsigned int currLoc=0; currLoc < simlocpos.size(); currLoc++)
    if(simlocpos[currLoc] > simloci)
      throw PrefExcept(LOCIPOS + " must have all positions within range of simulated number of loci " +
        itos(GetLociCount()));
  if(simlocpos.size() == 0 && num_simsets > 0)
    throw PrefExcept (LOCIPOS + " must be set for simulated data to be created");
  if(simAlleleLimits.size() > 0){
    if(simAlleleLimits.size() != 2)
      throw PrefExcept (SIMALLELELIMITS + " must have 2 frequencies listed");
    if(simAlleleLimits[0] < 0.0 || simAlleleLimits[1] > 0.5 || simAlleleLimits[1] < simAlleleLimits[0])
      throw PrefExcept (SIMALLELELIMITS + " must be between 0 and 0.5 and the first frequency must be smaller than the second");
  }
    /*
  if(simpopulation <= 0)
    throw PrefExcept(SIMPOPSIZE + " must be greater than 0");
  if(numgenerations <= 0)
    throw PrefExcept(SIMNUMGENS + " must be greater than 0");
  if(numgenes <= 0)
    throw PrefExcept(SIMGENES + " must be greater than 0");
  if(minnumsnp > maxnumsnp)
//std::cout << minnumsnp << " " << maxnumsnp << std::endl;
    throw PrefExcept(SIMMINSNP + " must be less than or equal to " + SIMMAXSNP);
  if(minnumsnp <= 0)
    throw PrefExcept(SIMMINSNP + " must be greater than 0");
  if(minrecomb > maxrecomb)
    throw PrefExcept(SIMMINRECOMB + " must be less than or equal to " + SIMMAXRECOMB);
  if(minrecomb <=0)
    throw PrefExcept(SIMMINRECOMB + " must be greater than 0");
	*/
}


string SimPreferences::itos(unsigned int number){
  stringstream oss;
  oss << number;
  return oss.str();
}

string SimPreferences::itos(long number){ stringstream oss;
  oss << number;
  return oss.str();
}

string SimPreferences::itos(int number){
  stringstream oss;
  oss << number;
  return oss.str();
}


string SimPreferences::itos(float number){
  stringstream oss;
  oss << number;
  return oss.str();
}

string SimPreferences::itos(double number){
  stringstream oss;
  oss << number;
  return oss.str();
}

// returns string with map element or throws KmodeExcept if doesn't
// exist
string SimPreferences::get_param(string key, map<string, string> & params){
        map<string, string>::iterator i = params.find(key);
        if(i == params.end()){
                string err = "Missing parameter " + key;
                throw PrefExcept(err);
        }
		//cout<<key<<":"<<i->second<<"\n";

        return i->second;  // returning string that is held in the map for the key
}


// skip rest of a line of input
void SimPreferences::skip_rest(ifstream & filestream){
        char c;
        do{
                c = filestream.get();
        }while(c != '\n' && c != EOF);
}

int SimPreferences::get_num_models()const{
   return num_models;
}

unsigned int SimPreferences::get_gen()const{
   return num_gens;
}

unsigned int SimPreferences::get_demes()const{
   return num_demes;
}

unsigned int SimPreferences::get_pop_size()const{
   return pop_size;
}

double SimPreferences::get_marg_var()const{
   return marg_var;
}

double SimPreferences::get_herit()const{
   return herit;
}

double SimPreferences::get_mutate()const{
   return mutation_rate;
}

double SimPreferences::get_cross()const{
   return crossover_rate;
}

double SimPreferences::get_table_var()const{
   return table_var;
}

int SimPreferences::get_update_interval()const{
   return update_interval;
}

bool SimPreferences::check_on_off(string status){
  if(status.compare("ON") == 0 || status.compare("on") == 0 || status.compare("On") == 0)
    return true;
  return false;
}

double SimPreferences::get_target_penetrance()const{
  return target_pen;
}

bool SimPreferences::datasimfile_set() const{
//cout << "simDataSimFile=" << simDataSimFile << endl;
  if(simDataSimFile.length() > 0)
    return true;
  else
    return false;
}

}
