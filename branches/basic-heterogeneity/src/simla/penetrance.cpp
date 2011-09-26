#include "penetrance.h"

namespace Simla {

Penetrance::Penetrance()
{
	Penetrance::maf = NULL;
	Penetrance::beta = NULL;
	Penetrance::models = NULL;
	Penetrance::size = 0;
	Penetrance::loci = 0;								// especially "loci" and "inter"(-actions) are important to be
	Penetrance::inter = 0;								// initialized to zero. The "startmodel()" function depends on it
	Penetrance::testindividuals = 0;
	Penetrance::prevalence = 0;
	Penetrance::rare = 1;										// rarest DX allele combination
	Penetrance::lociadded = 0;							// running count if we add loci via add_locus() function
	Penetrance::x0 = -999;
	Penetrance::beta = NULL;
	Penetrance::loci_list = NULL;
}

//***************************************

Penetrance::~Penetrance()
{
	delete [] models;
	delete [] beta;
	int i;
	for(i = 0; i < size; i++){
		delete [] loci_list[i];
	}
	delete [] loci_list;
	delete [] maf;
}

//***************************************
int Penetrance::init(const char *infile)
{
	int i;
	FILE *fptr = NULL;
	char str[100];
	char *tokenptr;
	fptr = fopen(infile, "r");
	if(fptr == NULL){
		perror("Error opening beta-file, exiting.");
		exit(0);
	}
	fgets(str, 100, fptr);
	fclose(fptr);
	
	Penetrance::loci = atoi(strtok(str, " "));
	Penetrance::inter = atoi(strtok(NULL, " "));
	Penetrance::maf = new double[loci];

	for(i = 0; i < loci; i++){
		maf[i] = atof(strtok(NULL, " "));
		rare *= maf[i];										// 2 alleles per locus
		rare *= maf[i];
	}
	if((rare <= 0) || (rare >= 1)){							// we require rare in (0,1)
		fprintf(stderr,"Penetrance::init. Problem with minor allele frequecies. Bailing out.\n");
		for(i = 0; i < loci; i++){
			fprintf(stderr,"%7.5f\n", maf[i]);
		}
		exit(0);
	}
	Penetrance::models = new double[loci];
	Penetrance::size = 0;
	
	for(i = 1; i <= inter; i++){
		size += n_over_k(loci, i);
	}
	Penetrance::beta = new double [size];
	Penetrance::loci_list = new int*[size];
	for(i = 0; i < size; i++){
		loci_list[i] = new int[loci];
		memset((void*)loci_list[i], -1, loci * sizeof(int));
	}
	read_data(infile);
	////// cdfbin ////////
	int which = 3;
	double p = certainty;
	double q = 1 - certainty;
	double s = min_rare;
	double xn;										// min number of individuals required
	double pr = rare;
	double ompr = 1 - pr;
	int status;
	double bound;
#ifdef WIN32
	cdfbin(&which, &q, &p, &s, &xn, &pr, &ompr, &status, &bound);
	if(status != 0){								// an error occured
		fprintf(stderr,"Penetrance::init. cdfbin() returned bad status of %i.\n", status);
        fprintf(stderr,"Penetrance::init. STATUS <-- 0 if calculation completed correctly\n               -I if input parameter number I is out of range\n                1 if answer appears to be lower than lowest search bound\n                2 if answer appears to be higher than greatest search bound\n                3 if P + Q .ne. 1\n                4 if PR + OMPR .ne. 1\n");
		fprintf(stderr,"Penetrance::init. p = %10.8f, q = %10.8f, s = %10.8f, xn = %12.2f, pr = %10.8f, ompr %10.8f\n", p, q, s, xn, pr, ompr);
		if(rare > 0){
			xn =  s / rare;
		}
		if((rare == 0) || (xn > max_test_individuals)){
			xn = max_test_individuals;
		}
	}
#else
		if(rare > 0){
			xn =  s / rare;
		}
		if((rare == 0) || (xn > max_test_individuals)){
			xn = max_test_individuals;
		}	
#endif
	Penetrance::testindividuals = (int)floor(xn+0.5);
	if(testindividuals > max_test_individuals){		// ... we don't want things to get out of hand in case there's some super rare allele combination
		testindividuals = max_test_individuals;
	}
	///// end cdfbin ///////
	return(testindividuals);
}

//***************************************
// ### input file format ###
// 5 3 0.1 0.1 0.1 0.1 0.1		// 5 disease loci, all possible 3-way interactions, DX allele freq. of 0.1 for each locus
// 0.025						// disease prevalence
// 0							// now come the disease models, 5 in this case
// 1							// 0 = recessive
// 0.5							// 0.5 = about multiplicative
// 0.5							// 1 = dominant
// 1
// 2 0.69						// locus 2 has a beta of 0.69
// 1-3-4 1.31					// interaction of loci 1,3,4 has beta of 1.31
// 5 0.69						// locus 5 has a beta of 1.31
// EOF							all others not listed are assumed to have a beta vaulue of zero

void Penetrance::read_data(const char *infile)
{
	FILE *fptr = NULL;
	fptr = fopen(infile, "r");
	int i,j;
	int index;
	int multiplier;
	int *loc = new int[loci];						// to store which loci are involved
	char str[1000];
	char* relevant;
	char *token;
	double bt;

	if(fptr == NULL){
		perror("Error opening beta-file, exiting.");
		exit(0);
	}
	fgets(str, 1000, fptr);							// skip "loci interaction"
	fgets(str, 1000, fptr);	
	token = strtok(str, " ");
	Penetrance::prevalence = atof(token);
	for(i = 0; i < loci; i++){						// first read disease modle for all loci
		fgets(str, 1000, fptr);
		token = strtok(str, " ");
		models[i] = atof(token);
	}
	memset((void*)beta, 0, size * sizeof(double));	// set all beta values to zero for now
	while(fgets(str, 1000, fptr)){
		memset((void*)loc, 0, loci * sizeof(int));	// reset all loci to zero
		relevant = strtok(str, " ");
		bt = atof(strtok(NULL, " "));				// the beta value
		if(bt != 0){								// allow negative beta values (protective "DX" locus)
			token = strtok(relevant, "-");
			while(token){
				j = atoi(token);
				if((j <= 0) || (j > loci)){
					fprintf(stderr,"Penetrance::read_data. Encountered invalid locus %i in control file. Bailing out\n", j);
					exit(0);
				}
				loc[j-1] = 1;
				token = strtok(NULL, "-");
			}
			multiplier = 1;
			index = get_index(loc);
			for(i = 0; i < loci; i++){
				loci_list[index][i] = loc[i];		// so we have a list of all loci involed in an interaction IF (beta > 0)
			}
			beta[index] = bt;
		}

	}
	
	fclose(fptr);
	delete [] loc;
}

//***************************************
// given loc[], like "0,0,1,0,1" we need to determine which index this corresponds too (here inteaction loci 3x5)
int Penetrance::get_index(int *loc)
{
	int i,j;
	int cnt = 0;
	int used = 0;
	int result = 0;
	int top = loci;
	int last = 1;
//int dbg_loc[6];
	int *list = new int[loci];
	for(i = 0; i < loci; i++){
//dbg_loc[i] = loc[i];
		if(loc[i] != 0){								// i.e == 1, but faster computational
			list[used] = i+1;
			used++;
		}
	}
	for(i = 1; i < used; i++){
		result += n_over_k(loci, i);
	}
	while(cnt < used){
		for(i = last, j = 1; i < list[cnt]; i++,j++){
			result += n_over_k(top-j, used-cnt-1);
		}
		last = list[cnt]+1;
		top = loci - list[cnt];
		cnt++;
	}

	delete [] list;
	return(result);

}

//***************************************

int Penetrance::n_over_k(int n, int k)
{
	if((k < 0) || (k > n)){
		return(0);
	}
	if((k == 0) || (n == k)){
		return(1);
	}
	int den = k;
	int num = n;
	int i;
	int result;
	for(i = 1; i < k; i++){
		num *= (n - i);
	}
	for(i = k-1; i > 1; i--){
		den *= i;
	}
	result = num / den;
	return(result);
}

//***************************************

void Penetrance::number2genotypes(int num, int *genotypes)
// [0] = genotype of locus 1
// ...
// [loci-1] = genotype of locus "loci"
// Example: num = 65
// 81	27	9	3	1
//  0	 2	1	0	2 = 2*27 + 1*9 + 0*3 + 2*1 = 65
// WE USE THE "BASE 3" SYSTEM 'COS EACH LOCUS CAN HAVE 0, 1 or 2 DX alleles
{
	memset((void*)genotypes, 0, loci * sizeof(int));					// set this mess to all zeros
	int multiplier = (int)powf(3.0, (float)(loci - 1));						//
	int remainder = num;
	int i = 0;
	if (num >= multiplier * 3){
		fprintf(stderr,"Penetrance::number2genotypes. Invalid genotype identifier of %i. Bailing out.\n", num);
		genotypes = NULL;
		i = *genotypes;
		exit(0);
	}
	while (remainder > 0 ){
		genotypes[i++] = remainder / multiplier;
		remainder = remainder % multiplier;		
		multiplier /= 3;
	}
}

//***************************************

double Penetrance::get_x0()
{
	return(Penetrance::x0);
}

//***************************************

double Penetrance::fx(int type)
{
	double arg = get_arg(type);
	arg += x0;
	double pen = exp(arg) / (1+exp(arg));
	return(pen);
}

//***************************************

double Penetrance::find_x0(int *genolist)
{
	double a;								// the x0-value to be...
	double upper = 10.0;
	double lower = -50.0;
	double current = -3.0;
	double epsilon = tolerance * prevalence;
	double accum;
	double temp;
	int i;
	int reps = 0;
	int done = 0;

	double *arglist = NULL;
	arglist = new double[testindividuals];
	for(i = 0; i < testindividuals; i++){
		arglist[i] = get_arg(genolist[i]);
	}
	while((done == 0) && (reps < 30)){
		accum = 0;
		for(i = 0; i < testindividuals; i++){
			temp = exp(current+arglist[i]) / (1+exp(current+arglist[i]));
			accum += temp;
		}
		accum /= (double)testindividuals;					// accum is now prevalence with "current" serving as x0
		if((accum >= (prevalence-epsilon)) && (accum <= (prevalence+epsilon))){
			a = current;									// within 1% tolerance relative to desired prevalence...Dat is good enough for gobernment work...
			done = 1;
		}
		else if(accum < prevalence){
			lower = current;
			current = (lower + upper) / 2.0;
		}
		else{												// obviously we have too high a prevalence 
			upper = current;
			current = (lower + upper) / 2.0;
		}
		reps++;
	}
	if(reps == 30){
		perror("x0 not converging\n");
		a = (lower + upper) / 2.0;
	}
	delete [] arglist;
	Penetrance::x0 = a;
	return(x0);
}

//***************************************

int Penetrance::get_loci()
{
	return(loci);
}

//***************************************

double Penetrance::get_arg(int type)
{
	double arg = 0;
	double temp;
	int i,j,k;
	int *indiv = new int[loci];
	double *indiv_model = new double[loci];
	number2genotypes(type, indiv);

	for(j = 0; j < loci; j++){					// need to transform the 0,1,2 genotypes into 0,dis_model,1 values used by logistic regression
		if(indiv[j] == 0){						// (1,1)
			indiv_model[j] = 0;
		}
		else if(indiv[j] == 2){					// (2,2)
			indiv_model[j] = 1;
		}
		else{									// (1,2) == (2,1) heterozygot case
			indiv_model[j] = models[j];			// disease model for locus as specified in control file
		}
	}

	for(j = 0; j < size; j++){					// over all beta values
		if(beta[j] != 0){						// otherwise we can care less
			temp = beta[j];
			for(k = 0; k < loci; k++){
				if(loci_list[j][k] == 1){		// if the kth locus is part of this interaction...proceed
					temp *= indiv_model[k];		// indiv_model[] releates to ONE individual's genotype at the k-disease loci. 0=(1,1), 1=(2,2), model=(1,2) genotype
				}
			}
			arg += temp;
		}
	}
	delete [] indiv;
	delete [] indiv_model;
	return(arg);
}

//***************************************

void Penetrance::start_model(int loc, int interactions, double prev)
{
	int i;
	if((loc < 1) || (loc > max_loci)){
		fprintf(stderr,"Penetrance::start_model. Invalid number of DX loci (%i). Try again.\n", loc);
		return;
	}
	if((interactions < 1) || (interactions > max_interactions)){
		fprintf(stderr,"Penetrance::start_model. Invalid number of interactions (%i). Try again.\n", interactions);
		return;
	}
	if((prev <= 0) || (prev >= max_prevalence)){
		fprintf(stderr,"Penetrance::start_model. Invalid disease prevalence of (%7.5f). Try again.\n", prev);
		return;
	}
	if(Penetrance::loci == 0){
		Penetrance::loci = loc;
		Penetrance::prevalence = prev;
		Penetrance::inter = interactions;
		Penetrance::maf = new double[loci];
		Penetrance::models = new double[loci];
		for(i = 1; i <= inter; i++){
			size += n_over_k(loci, i);
		}
		Penetrance::beta = new double [size];
		memset((void*)beta, 0, size * sizeof(double));			// so all beta values are 0, unless you tell us otherwise in "add_interaction()"
		Penetrance::loci_list = new int*[size];
		for(i = 0; i < size; i++){
			loci_list[i] = new int[loci];
			memset((void*)loci_list[i], -1, loci * sizeof(int));
		}
	}
	else{
		fprintf(stderr,"Penetrance::start_model. This function can only be called once. Bailing out\n");
		exit(0);
	}
}

//***************************************

void Penetrance::add_locus(double model, double freq)
{
	if(lociadded < loci){
		if((freq <= 0) || (freq >= 1)){
			fprintf(stderr,"Penetrance::add_locus. Invalid allele frequency of %7.5f. Bailing out.\n", freq);
			exit(0);
		}
		if((model < 0) || (model > 1)){
			fprintf(stderr,"Penetrance::add_locus. Invalid disease model of %7.5f. Bailing out.\n", freq);
			exit(0);
		}
		maf[lociadded] = freq;
		rare *= freq;
		rare *= freq;										// 2 alleles per person
		models[lociadded] = model;
		lociadded++;
	}
	else{
		fprintf(stderr,"Penetrance::add_locus. Attempt to add past the maximum number of loci ignored.\n");
	}
}

//***************************************

void Penetrance::add_interaction(int members, int *inter, double bt)
{
	if((members < 1) || (members > Penetrance::inter)){
		fprintf(stderr,"Penetrance::add_interaction. An interaction of %i loci is not possible. Bailing out.\n", members);
		exit(0);
	}
	int *loc = new int[loci];
	int i;
	int index;
	memset((void*)loc, 0, loci * sizeof(int));
	for(i = 0; i < members; i++){
		loc[inter[i]-1] = 1;
	}
	index = get_index(loc);
	for(i = 0; i < loci; i++){
		loci_list[index][i] = loc[i];
	}
	beta[index] = bt;
	delete [] loc;
}

//***************************************

int Penetrance::close_model()
{
	if(lociadded != loci){
		fprintf(stderr, "Penetrance::close_model. Loci added = %i, total loci = %i. Can not close until all loci are added.\n", lociadded, loci);
		return(0);
	}
	////// cdfbin ////////
	int which = 3;
	double p = certainty;
	double q = 1 - certainty;
	double s = min_rare;
	double xn;										// min number of individuals required
	double pr = rare;
	double ompr = 1 - pr;
	int status;
	double bound;
//I'm getting rid of this, since it isn't working under gcc 3.2 (mingw) and wasn't in use for linux
#ifdef __THIS_IS_NOT_DEFINED__// WIN32
	cdfbin(&which, &q, &p, &s, &xn, &pr, &ompr, &status, &bound);
	if(status != 0){								// an error occured
		fprintf(stderr,"Penetrance::close_model. cdfbin() returned bad status of %i.\n", status);
        fprintf(stderr,"Penetrance::close_model. STATUS <-- 0 if calculation completed correctly\n               -I if input parameter number I is out of range\n                1 if answer appears to be lower than lowest search bound\n                2 if answer appears to be higher than greatest search bound\n                3 if P + Q .ne. 1\n                4 if PR + OMPR .ne. 1\n");
		fprintf(stderr,"Penetrance::close_model. p = %10.8f, q = %10.8f, s = %10.8f, xn = %12.2f, pr = %10.8f, ompr %10.8f\n", p, q, s, xn, pr, ompr);
		if(rare > 0){
			xn =  s / rare;
		}
		if((rare == 0) || (xn > max_test_individuals)){
			xn = max_test_individuals;
		}
	}
#else
		if(rare > 0){
			xn =  s / rare;
		}
		if((rare == 0) || (xn > max_test_individuals)){
			xn = max_test_individuals;
		}	
#endif
	Penetrance::testindividuals = (int)floor(xn+0.5);
	if(testindividuals > max_test_individuals){		// ... we don't want things to get out of hand in case there's some super rare allele combination
		testindividuals = max_test_individuals;
	}
	else if(testindividuals < min_test_individuals){
		testindividuals = min_test_individuals;
	}
	///// end cdfbin ///////
	return(testindividuals);
}

//***************************************

void Penetrance::reset()
{
	delete [] models;
	delete [] beta;
	int i;
	for(i = 0; i < size; i++){
		delete [] loci_list[i];
	}
	delete [] loci_list;
	delete [] maf;
	Penetrance::maf = NULL;
	Penetrance::beta = NULL;
	Penetrance::models = NULL;
	Penetrance::size = 0;
	Penetrance::loci = 0;								// especially "loci" and "inter"(-actions) are important to be
	Penetrance::inter = 0;								// initialized to zero. The "startmodel()" function depends on it
	Penetrance::testindividuals = 0;
	Penetrance::prevalence = 0;
	Penetrance::rare = 1;										// rarest DX allele combination
	Penetrance::lociadded = 0;							// running count if we add loci via add_locus() function
	Penetrance::x0 = -999;
	Penetrance::beta = NULL;
	Penetrance::loci_list = NULL;
}

//***************************************

}//Simla

//***************************************
//***************************************
