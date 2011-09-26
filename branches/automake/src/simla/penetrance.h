#ifndef __PENETRANCE_H__
#define __PENETRANCE_H__

#include <math.h>
#include <stdlib.h>
#include "random.h"
#include "cdflib.h"


namespace Simla {

const double tolerance = 0.001;
const double certainty = 0.95;								// 95% certainty to have at least
const double min_rare = 5.0;								// "min_rare" given all DX loci allele frequencies (i.e. the rarest allele combination over all DX loci).
const double max_prevalence = 0.3;							// maximum disease prevalence. Be careful, high prevalences lead to unexpected results. In the extreme consider a prevalence of 1. Can you still have a relative risk > 1??? ...no
const int max_test_individuals = 5000000;					// the most "random individuals" we can ask for when determining x0
const int min_test_individuals = 5000;						// absolute minimum nuber of individuals to test
const int max_loci = 6;										// max number of DX (disease) loci
const int max_interactions = 6;								// max number of interactions for now...probably can increase this number to about 6 or 7 (but <= loci) without causing a crash

// Quick tutorial.
// For this to work you need a control file that looks as follows:
/* ### EXAMPLE CONTROL FILE ###
6 3	0.1 0.2 0.3 0.4 0.5 0.6			// 6 disease loci and up to 3-way intereaction. Respective DX allele frequencies are: 0.1, 0.2, ...
0.05								// targeted disease prevalence of 0.05
0									// the next "dx loci, in this example 6" entires give the
0.5									// disease model for each (of the 6 in this example) disease locus.
1									// 0 = recessive, 1 = dominant, 0.5 ~ multiplicative
1									// any real number in [0,1] is valid, the higher, the more "dominant"
0									// the lower the more "recessive" the model becomes.
0.5
2        0.02						// loci combinations followed associted BETA value					
4        0.04						// Any combination not listed is assumed to have a beat value of 0.
1-5      0.1						// The values in this example were useful for debugging but are probably
1-2-5    0.24						// pretty meaningless when it comes to investigation genetic models.
1-3-6    0.28
2-3-4    0.32
3-4-6    0.39
4-5-6    0.41

Important notes: 
1. DO NOT have any line starting with whitespace. The reader expects the first character to be 0-9.
2. Fields are seperated by one or more whitespaces.
3. The loci interactions/combinations are one continous string with "-" seperating the individual loci.
4. The loci combinations don't have to be in any particular order. i.e. it doesn't matter if 1-2-5 comes before or after 2-3-4
5. Within one combination the loci must be ordered. i.e. 2-4-5 is the only acceptable format of this 3-loci interaction. 4-2-5 for example is not acceptable.
6. The last line of the control file is a newline by itself.
*/
//	to get things started:
// 1. Declare "Penetrance-Class" variable.
// 2. Call Penetrance::init("name-of-your-control-file-here")
// 3. It will call an external function "generate_genotypes(int number_of_loci, int test_sample_size, int *test_list)"
//	- the memory for test_list has been allocated. The generate_genotypes() function is supposed to 
//	fill in the list with "test_sample_size" RANDOM genotypes like they are likely to occur in the
//	general population. See "void generate_genotypes(int num_of_loci, int n, int *list)" as example. Each bi-allelic 
//	locus (the only type we consider) has 3 possible genotypes. So with n-loci we have a total of 3^n
//	distinct disease genotypes, coded as integers 0 .. (3^n - 1).
// 4. To get penetrance values call; double Penetrance::fx(int individuals-genotype)
//	The argument "individuals-genotype" is an integer coded as in the generate_genotypes() function.
//	The return value is the penetrance value in (0,1). Usually you would compare this value
//	to a unifom random number on [0,1] to decide affection status. If random less than penetrance => affected
//	else unaffected.
class Penetrance{

public:
	Penetrance();
	~Penetrance();
	int init(const char *infile);							// returns number of random individuals to "find_x0()"		
	int n_over_k(int n, int k);						// exactly what you think it is
	void read_data(const char *infile);					// ...got to get the data from somewhere...
	int get_index(int *loc);						// inverse of above, given the loci involved, tell me the index. 
	double find_x0(int*);							// the all important X0, needed for the penetrance function
	void number2genotypes(int num, int *genotypes); // give me a number and I get you the genotypes, 0 = 1,1; 1 = 1,2 or 2,1; 2 = 2,2. Always based on "loci" total loci.
	double get_arg(int type);						// put in number coded genotype, returns arg for individual
	double fx(int type);							// arg: the integer coded genotype of an individual, returns penetrance value.
													// It is up the the calling function to interpret the penetrance value.
	int get_loci();									// returns number of DX loci
	double get_x0();								// returns x0
	// 
	///// the following functions can be used INSTEAD of the init() function which takes a control file ////
	// The calling sequence is something like this:
	// start_model(int number_of_disease_loci, int number_of_loci_in_largest_interaction, double prevalence) ... which is at minimum 1 and at maximum some very finite number >= 3. So I only test it with 1, 2 and 3. When you get a segfault you know ...
	// add_locus(double model, double DX_allele_frequency);
	// ..., call add_locus exactly once for each locus specified via start_model();
	// add_interaction(int members, int *with_id_of_each_member_locus_of_interaction);  locus numbering is 1,2,...,loci
	// ..., call as many times as you like. (members >= 1) && (members <= number_of_loci_in_largest_interaction)
	// NOTE!!! Here you MUST also add the beta values of single loci. Think of them as interacting with themselves or whatever
	// int close_model();			finalizes the setup. returns number of needed random test specimens to determine x0.
	void start_model(int loci, int interactions, double prevalence);	
	void add_locus(double model, double freq);
	void add_interaction(int members, int *inter, double beta);
	int close_model();
	void reset();									// In case you want to start over again

private:
	double x0;										// Yes, THE x0 used in logistic regression function
	double prevalence;								// disease prevalence
	int loci;										// number of disease loci
	int lociadded;									// running count used by ad_locus() function
	int inter;										// up to & including interactions of degree "inter"
	int size;										// number of beta values needed to cover all main effects + interactions
	double rare;									// rarest allele combination. Relevant for how many test individuals we need
	double *beta;									// beta-values
	double *models;									// disease models for loci
	int testindividuals;							// number returned by "int init()". "double find_x0" expects that many random individuals
	double *maf;									// "minor_allele_frequencies". i.e. allele freq. for DX allele at each DX locus
	int **loci_list;								
	// lists the loci involved for the RELEVANT (i.e. beta != 0) interactions
	// [0] = -1,...,-1  with -1 indication that beta for this interaction equals zero nad we do not have any info which loci are involved (even so it is obvious in this case that only locus 1 would be involved)
	//  ...
	// [k] = 0,1,0,1,0,0 in this 6-locus model we have loci 2 and 4 interacting as indicated by the 1s in the appropriate places
	// This info is important when we calculate thepenetrance values.
	// if (beta == 0)...fine, nothing to do...NEXT! Note: beta >= 0...always. We don't deal with negative beta values.
	// if (beta > 0) now we need to know which loci are involved, so we look at "0,1,0,1,0,0" and know that this beta value
	// corresponds to loci 2x4. So, according to the individuals genotypes & disease models @ loci 2 & 4 we can calculate
	// the arg for the logisitc regression function.
};

} // end namespace SIMLA
#endif
