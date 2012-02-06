#include "main.h"
#include "externals.h"
using namespace Simla;

int main(int argc, char *argv[]){
	char infile[100];
	if(argc >= 2){
		strcpy(infile, argv[1]);
	}
	Random::setseed(1371);
	Random r;
	r.init();
	Penetrance p;
	double x0;
//// Now let's test stuff /////
	int i;
	
							// setup data via function calls
	p.start_model(3, 2, 0.25);			// start model with args: loci, max_size_interaction, DX_prevalence
	cout<<"Initializing 3 loci with disease at 0.48 frequency\n";
	p.add_locus(1.0, 0.48);				// add locus with "DX model [0,1]" and DX allele frequency 
	p.add_locus(0.5, 0.48);				// args: DX_model in [0,1], DX_allele_frequency
	p.add_locus(0.5, 0.48);
	int *llist = new int[2];			// allocate enough to take care of largest interaction....here 2
	llist[0] = 1;						// enter the locus/loci involved					
	p.add_interaction(1, llist, 0.0);	// args: number_of_loci_in_interaction, array_containing_loci_ID(s), beta
	llist[0] = 2;
	p.add_interaction(1, llist, 0.0);
	llist[0] = 3;
	p.add_interaction(1, llist, 0.0);
	llist[0] = 2;
	llist[1] = 3;
	p.add_interaction(2, llist, 1.6);
	int samples = p.close_model();		// close model and get number of samples required to find x0

	int *random_folks = new int[samples];
	double *freq = new double[3];		// 3 == number of DX-loci
	freq[0] = 0.48;
	freq[1] = 0.48;
	freq[2] = 0.48;
	
	generate_genotypes(3, freq, samples, random_folks);	// this batch is to find x0
	p.find_x0(random_folks);
	delete [] random_folks;

	random_folks = new int[100000];
	generate_genotypes(3, freq, 100000, random_folks);	// the 2nd batch to test how "good" an x0 we found
	delete [] freq;



	int tally = 0;
	double penetr;
	double rand;
	for(i = 0; i < 100000; i++){
		penetr = p.fx(random_folks[i]);
		rand = r.ran1();
		if(rand < penetr){
			tally++;
		}
	}
	delete [] random_folks;
	penetr = (double)tally / 100000.0;
	printf("observed disease prevalence is %8.6f\n",penetr);
	x0 = p.get_x0();
	printf("Based on %i test samples and an x0 value of %8.6f\n",samples, x0);
////////// print out penetrances of all combinantions //////////////////
	int allcombinations = 1;
	int loci = p.get_loci();
	int j;
	double penetrance;
	int *indiv = new int[loci];

	for(i = 0; i < loci; i++){
		allcombinations *= 3;					// 3 genotypes per locus
	}
	for(i = 0; i < allcombinations; i++){
		penetrance = p.fx(i);
		p.number2genotypes(i, indiv);
		printf("%i",indiv[0]);
		for(j = 1; j < loci; j++){
			printf("-%i",indiv[j]);
		}
		printf("   %7.5f\n",penetrance);
	}
	
	///// test reset //////////
	// 2nd test///////////
	
	r.init();
	p.reset();
							// setup data via function calls
	p.start_model(3, 2, 0.25);			// start model with args: loci, max_size_interaction, DX_prevalence
	cout<<"initializing 3 loci with disease at 0.52 frequency\n";
	p.add_locus(1.0, 0.52);				// add locus with "DX model [0,1]" and DX allele frequency 
	p.add_locus(0.5, 0.52);				// args: DX_model in [0,1], DX_allele_frequency
	p.add_locus(0.5, 0.52);

	llist[0] = 1;						// enter the locus/loci involved					
	p.add_interaction(1, llist, 0.0);	// args: number_of_loci_in_interaction, array_containing_loci_ID(s), beta
	llist[0] = 2;
	p.add_interaction(1, llist, 0.0);
	llist[0] = 3;
	p.add_interaction(1, llist, 0.0);
	llist[0] = 2;
	llist[1] = 3;
	p.add_interaction(2, llist, 1.6);
	samples = p.close_model();		// close model and get number of samples required to find x0

	random_folks = new int[samples];
	freq = new double[3];		// 3 == number of DX-loci
	freq[0] = 0.48;
	freq[1] = 0.48;
	freq[2] = 0.48;

	generate_genotypes(3, freq, samples, random_folks);	// this batch is to find x0
	p.find_x0(random_folks);
	delete [] random_folks;
	random_folks = new int[100000];
	generate_genotypes(3, freq, 100000, random_folks);	// the 2nd batch to test how "good" an x0 we found
	delete [] freq;
/////////////// now test the 2nd model ///////////
	tally = 0;
		for(i = 0; i < 100000; i++){
		penetr = p.fx(random_folks[i]);
		rand = r.ran1();
		if(rand < penetr){
			tally++;
		}
	}
	delete [] random_folks;
	penetr = (double)tally / 100000.0;
	printf("observed disease prevalence is %8.6f\n",penetr);
	x0 = p.get_x0();
	printf("Based on %i test samples and an x0 value of %8.6f\n",samples, x0);
/////////////////////////////
	allcombinations = 1;
	loci = p.get_loci();

	indiv = new int[loci];
	for(i = 0; i < loci; i++){
		allcombinations *= 3;					// 3 genotypes per locus
	}
	for(i = 0; i < allcombinations; i++){
		penetrance = p.fx(i);
		p.number2genotypes(i, indiv);
		printf("%i",indiv[0]);
		for(j = 1; j < loci; j++){
			printf("-%i",indiv[j]);
		}
		printf("   %7.5f\n",penetrance);
	}
/////////////////////////////
	delete [] indiv;
	delete [] llist;
////////////////////////////
	return(0);
}
