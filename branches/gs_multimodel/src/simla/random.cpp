#include "random.h"

#if defined (WIN32)
	#include <windows.h>
//	#include <afxwin.h>
#endif

namespace Simla {

static time_t testvar = time(NULL);
static long idum = (testvar < 0) ? (long) testvar : (long) -1*testvar;
static long s1 = idum;
static long s2 = 7701;
uint32 init_key[def_init] = {0,0,0,0};
static int inits = 0;
static uint32 seed = 0;
// Knuth related
double ranf_arr_buf[QUALITY];
double ranf_arr_dummy=-1.0, ranf_arr_started=-1.0;
double *ranf_arr_ptr=&ranf_arr_dummy; /* the next random fraction, or -1 */
double ran_u[KK];           /* the generator state */

// end Knuth related

int Random::maxbits = 0;
uint32 Random::maxint = 1;
int Random::negmaxbits = 0;
uint32 Random::mt[NN] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 
	}; 
static int mti=NN+1; /* mti==N+1 means mt[N] is not initialized */

//***************************************************

Random::Random()
{
    
}

//****************************************************

void Random::init()
{
  static long repet = 0;
  static uint32 agitator = 0xb205a219;

  uint32 utemp;
  repet++;
  static int ini = 0;
  ifstream inf;
  inf.open("randgs.txt",ios::in);
 
  uint32 a = 1;
  uint32 b = -1;							// b = 1...1 regardless of number of bits
  b >>= 1;								// b = 01...1
  b = ~b;								// b = 10...0
  Random::maxbits = 0;
  maxint = 1;
  for(; a !=0; a <<= 1, b >>= 1){
	utemp = maxint;
	maxint <<=1;
	maxint |= utemp;	
	maxbits++;
  }
  if(maxbits != 32){
	cerr <<"Random::init. 32-bit unsigned data type has "<<maxbits<<"-bits. Bailing out.\n";
	exit(0);
  }
  char line[200];
  int i;
  char *tokenptr;

/*
		char *tokenptr;
		
		infile.getline(line, max_rdln);
		tokenptr = strtok(line, " ");
		ptr_to_buffer(tokenptr, name, max_length);
		tokenptr = strtok(NULL, " ");
		type = *tokenptr;                                 // 'X' or 'A'
		tokenptr = strtok(NULL, " ");
		alleles = atoi(tokenptr);
		tokenptr = strtok(NULL, " ");
		loc = atof(tokenptr);
		tokenptr = strtok(NULL, " ");
  */
//cerr <<"using init()"<<endl;
  if(inf.is_open()){
	  inf.getline(line, 200);
	  tokenptr = strtok(line, " ");
	  i = 0;
	while((tokenptr != NULL) && (i < def_init)){
		utemp = (unsigned long)atoi(tokenptr);
		utemp = shuffle_seed(utemp);
		init_key[i] = utemp;
		tokenptr = strtok(NULL, " ");
		i++;
	}
	inf.close();
	if(i == 4) ini++;
  }
  if(!ini){
  
  
#if defined (OSX) 
	uid_t u;
	u = getuid();

	long sum = (long)u;
	sum += (long)(getsid((pid_t)0));
	testvar = time(NULL);
  	long temp = (long)testvar;
	temp %= 26617;
	long pid = (long)(getpid() % 130817);
	idum = (unsigned int)((pid * temp + sum + temp + pid ) % 9217433);
	s2 = (s2 + s2 * temp + repet + idum) % 2502133;
	seed = (uint32)s2;
	seed = shuffle_seed(seed);
	seed ^= agitator;
	seed = shuffle_seed(seed);

#endif

#if defined (SOLARIS)  
	uid_t u;
	u = getuid();
	hrtime_t h;
	h = gethrtime();
	long hr = (long)h;
	hr %= 659411;
	if(hr < 0){
		hr *= -1;
	}
	long sum = (long)u;
	sum += (long)(getsid((pid_t)0));
	testvar = time(NULL);
  	long temp = (long)testvar;
	temp %= 26617;
	long pid = (long)(getpid() % 130817);
	idum = (unsigned int)((pid * temp + sum + temp + pid + hr) % 9217433);
	s2 = (s2 + s2 * temp + repet + idum) % 2502133;	
	seed = (uint32)s2;
	seed = shuffle_seed(seed);
	seed ^= agitator;
	seed = shuffle_seed(seed);
	struct timeb tp;
	uint32 millisec;
	if(0 == ftime(&tp)){
	  millisec = (uint32)tp.millitm;
	  millisec *= 17;
	  seed ^= millisec;
	  seed = shuffle_seed(seed);
	}

#endif

#if defined (LINUX)
	uid_t u;
	u = getuid();
	int h;
	struct rusage tmp;
	long hr = 0;
	long accum = 0;
#ifdef RUSAGE_SELF                // this introduces some more randomness into the equation
	getrusage(RUSAGE_SELF, &tmp);
	hr += (long)tmp.ru_utime.tv_sec;
	hr += (long)tmp.ru_utime.tv_usec;
	hr += tmp.ru_maxrss;
	hr += tmp.ru_ixrss;
	hr += tmp.ru_idrss;
	hr += tmp.ru_minflt;
	hr += tmp.ru_majflt;
	hr += tmp.ru_nswap;
	hr += tmp.ru_isrss;
	hr += tmp.ru_inblock;
	getrusage(RUSAGE_SELF, &tmp);
	accum += (long)tmp.ru_utime.tv_sec;
	accum += (long)tmp.ru_utime.tv_usec;
	accum += tmp.ru_maxrss;
	accum += tmp.ru_ixrss;
	accum += tmp.ru_idrss;
	accum += tmp.ru_minflt;
	accum += tmp.ru_majflt;
	accum += tmp.ru_nswap;
	accum += tmp.ru_isrss;
	accum += tmp.ru_inblock;
	hr += (accum / 1000) * 17;
	hr += (accum / 100) * 13;
	hr += (accum / 10) * 7;
	hr += accum % 23;

#else
	hr = -1;
#endif

	hr %= 659411;
	if(hr < 0){
		hr *= -1;
	}
	long sum = (long)u;
	sum += (long)(getsid((pid_t)0));
	testvar = time(NULL);
  	long temp = (long)testvar;
	temp %= 26617;
	long pid = (long)(getpid() % 13081);
	idum = (unsigned int)((pid * temp + sum + temp + pid + hr) % 9217433);
	s2 = (s2 + s2 * temp + repet + idum) % 2502133;
	seed = (uint32)s2;
	seed = shuffle_seed(seed);
	seed ^= agitator;
	seed = shuffle_seed(seed);

#endif

#if defined (WIN32)
	
	testvar = time(NULL);
	uint32 temp = (uint32)testvar;
	SYSTEMTIME st;                  // Can't find appropriate header file so that this will compile.
//	CTime time;
	GetSystemTime(&st);
	uint32 millisec = (uint32)(17 * st.wMilliseconds);

		testvar = time(NULL);
	 temp %= 26617;
	 long pid = (long)(_getpid() % 13081);
	 srand(temp + pid);
	 long r = rand();
	 long s = (long)(&r) + (long)(&init) + (long)(&pid);    // we take randomness wherever we can find it...
	 idum = (unsigned int)((pid * temp + temp + pid + r + s) % 9217433);
	 seed = (uint32)(s2 + idum + repet + s2 * temp);
	seed = shuffle_seed(seed);
	seed ^= millisec;
	seed = shuffle_seed(seed);
	seed ^= agitator;
	seed = shuffle_seed(seed);


	/*
	uint32 pid = (uint32)_getpid();
	 long longseed = (long)(temp + millisec);
	 srand(longseed);
	 uint32 r = rand();

	 seed = shuffle_seed(temp);
	 seed ^= millisec;
	 seed = shuffle_seed(seed);
	 seed ^= r;
	 seed = shuffle_seed(seed);
     seed ^= pid;
	 seed = shuffle_seed(seed);
	 seed ^= (uint32)repet;
	 seed = shuffle_seed(seed);
	 seed ^= agitator;
	 seed = shuffle_seed(seed);
*/

#endif
	ini = 1;
	mti = NN + 1;
	init_key[0] = (unsigned long)seed;
	init_key[1] = (unsigned long)idum;
	utemp = idum >> 3;
	utemp = seed ^ utemp;
	init_key[2] = (unsigned long)utemp;
	utemp = seed << 4;
	utemp = idum ^ utemp;
	utemp = shuffle_seed(utemp);
	init_key[3] = (unsigned long)utemp;
  }
  init_by_array(init_key);
  for(i = 0; i < def_init; i++){
    seed ^= init_key[i];
  }
  ranf_start(seed);
}

//***************************************************

Random::~Random()
{

}

//***************************************************

void Random::setseed(uint32 i)
{
  int j;
  uint32 utemp;
  seed = i;
  bitsinunsignedlong();
  if((seed == 0) && (inits == 0)){
    //cerr <<"seed was set to zero"<<endl;
    init();
    inits = 1;
  }
  else if(inits == 0){
    utemp = i;
    for(j = 0; j < def_init; j++){
      utemp &= 0xfff00fff;
      utemp = shuffle_seed(utemp);
      utemp |= 0x000ff000;
      utemp = shuffle_seed(utemp);
	  init_key[j] = utemp;
	  utemp++;
    }
    init_by_array(init_key);
    mti = NN + 1;
    inits = 1;
    // for Knuth's generator
  for(i = 0; i < (int)def_init; i++){
    seed ^= init_key[i];
  }
  ranf_start(seed);
    //cerr <<"used setseed, seed = "<<i<<endl;
  }
}

//***************************************************

void Random::bitsinunsignedlong()
{
	maxbits = 0;
  uint32 max = 0;
  for(max--; max != 0; maxbits++) max >>= 1;	// find number of bits in type "unsigned long"
  negmaxbits = -1 * maxbits;

}

//***************************************************

uint32 Random::testseed(){
	seed = 0;
	inits = 0;
	setseed(0);
	return(seed);
}

//***************************************************

void Random::clearseed()
{
	seed = 0;
	inits = 0;
}

//***************************************************

int Random::ranint(int arg)
{
	int result;
	uint32 rand = ran2();
//	rand >>= 3;						// to eliminate possible patterns in the data
	result = (int)(rand % arg);
	return(result);
}

//****************************************************

double Random::rnorm(double mu, double dev)
{
	
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	double result;
	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		result = v2*fac;
	} else {
		iset=0;
		result = gset;
	}
	result *= dev;
	result += mu;
	return(result);
}	

//***************************************************

double Random::ranb()
{
// This library is free software; you can redistribute it and/or   
// modify it under the terms of the GNU Library General Public     
// License as published by the Free Software Foundation; either    
// version 2 of the License, or (at your option) any later         
// version.                                                        
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of  
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            
// See the GNU Library General Public License for more details.    
// You should have received a copy of the GNU Library General      
// Public License along with this library; if not, write to the    
// Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA    
// 02111-1307  USA                                                 

// Copyright (C) 1997, 1999, 2001                                  
//    Makoto Matsumoto and  Takuji Nishimura.                      
//                                                                 
// Any feedback is very welcome. For any question, comments,       
// see http://www.math.keio.ac.jp/matumoto/emt.html or email       
// matumoto@math.keio.ac.jp                                        

// REFERENCE                                                       
// M. Matsumoto and T. Nishimura,                                  
// "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  
// Pseudo-Random Number Generator",                                
// ACM Transactions on Modeling and Computer Simulation,           
// Vol. 8, No. 1, January 1998, pp 3--30.            
	
// The quality of this generator is very good. Unfortunately is is not the fastest.
// On a particular test with low disease prevalence SIMLA slowed down from
// 3:20 to 13:50 minutes when compared to either of Knuth's generators.
	double result;
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= NN) { /* generate N words at one time */
        int kk;

        if (mti == NN+1)   /* if sgenrand() has not been called, */
            sgenrand(seed); /* a default initial seed is used   */

        for (kk=0;kk<NN-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);					// y >> 11
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;		// y << 7  & 0x9d2c5680
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;		// y << 15 & 0xefc60000
    y ^= TEMPERING_SHIFT_L(y);					// y >> 18
	result = (double)y;
	result = ldexp(result, negmaxbits);
	return(result);
//    return ( ((double)y + 1.0) * 2.3283064359965952e-10 ); /* reals: (0,1)-interval */
    /* return y; */ /* for integer generation */
}

//***************************************************

uint32 Random::ran2()
{
     return((uint32)floor(ran1() * maxint));

}

//***************************************************

double Random::ranc()
//Knuth, D.E. 1981 Seminumerical Algorithms, 2nd ed. vol.2 of The Art of Computer Programming
// (Reading, MA: Addison-Wesley), 3.2-3.3. [4]
// The quality of this generator is decent but not great. It is extremely fast.
{
	static int inext, inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if(iff == 0){
		iff = 1;
		mj=labs(MSEED-labs(seed));
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for(i = 1; i <= 54; i++){
			ii = (21*i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if(mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for(k = 1; k <= 4; k++){
			for(i = 1; i <= 55; i++){
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		}
		inext = 0;
		inextp = 31;
//		seed = 1;                    // indicates that we are done with initialization
	}
	// Here is where the real algorithm starts
	if (++inext == 56) inext = 1;
	if (++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return (mj*FAC);
//	return result;
}

//***************************************************

void Random::bye()
{
	ofstream out;
	uint32 temp;
	uint32 temp1;
	uint32 temp2;
	int i;
	out.open("randgs.txt", ios::out);
	for(i = 0; i < 4; i++){
		temp1 = ran2();
		temp1 = init_key[i] ^ temp1;
		temp1 = shuffle_seed(temp1);
		temp2 = ran2();
		temp2 = init_key[i] ^ temp1;
		temp2 = shuffle_seed(temp1);
		temp = temp1 ^ temp2;
		out <<temp<<"  ";
	}
	out <<endl;
	out.close();
}

//***************************************************

double Random::normdist(double mu, double dev, double x)
{
	double result = 0;
	double temp1 = 0;
	double coeff;
	double sqr;
	double expo;
	temp1 = dev * sqrt2pi;
	coeff = 1 / temp1;
	sqr = (x-mu)*(x-mu);
	expo = -1/(2*dev*dev);
	expo *= sqr;
	result = coeff * exp(expo);
	return(result);
}

//***************************************************

void Random::sgenrand(unsigned long seed)
{
    int i;
	
    for (i=0;i<NN;i++) {
         mt[i] = seed & 0xffff0000;
         seed = 69069 * seed + 1;
         mt[i] |= (seed & 0xffff0000) >> 16;
         seed = (69069 * seed + 1) % 82376621; // is a prime number
    }
}

//***************************************************

uint32 Random::shuffle_seed(uint32 seed)
{
	uint32 hash = seed ^ 0x96c3a5d2;
	uint32 temp;
	uint32 a = 1;
	uint32 b;
	uint32 reverse = 0;
	uint32 result = 0;

	uint32 first;
	uint32 firstb;

	uint32 second;
	uint32 secondb;

	uint32 third;
	uint32 thirdb;

	uint32 fourth;
	uint32 fourthb;

	uint32 fifth;
	uint32 fifthb;

	uint32 sixth;
	uint32 sixthb;

	uint32 seventh;
	uint32 seventhb;

	uint32 eighth;
	uint32 eighthb;

// reverse the bits in "hash"
// init a = 0...01
// init b = 10...0 problem is we dont know how many bits we got 16, 32, 64?
	b = -1;							// b = 1...1 regardless of number of bits
	b >>= 1;						// b = 01...1
	b = ~b;							// b = 10...0

	for(; a !=0; a <<= 1, b >>= 1){
		temp = b & hash;
		if(temp){
			reverse |= a;
		}
	}

		first = reverse & 0xc0000000;
		firstb = reverse & 0x30000000;

		second = reverse & 0x0c000000;
		secondb = reverse & 0x03000000;

		third = reverse & 0x00c00000;
		thirdb = reverse & 0x00300000;

		fourth = reverse & 0x000c0000;
		fourthb = reverse & 0x00030000;

		fifth = reverse & 0x0000c000;
		fifthb = reverse & 0x00003000;

		sixth = reverse & 0x00000c00;
		sixthb = reverse & 0x00000300;

		seventh = reverse & 0x000000c0;
		seventhb = reverse & 0x00000030;

		eighth = reverse & 0x0000000c;
		eighthb = reverse & 0x00000003;


	// now put it back together as 4b-7-8b-2-3b-6b-1-3-5b-8-1b-7b-2b-6-4-5
	first >>= 12;
	second >>= 2;
	third >>= 6;
	fourth >>= 16;
	fifth >>= 14;
	sixth >>= 6;
	seventh <<= 22;
	eighth <<= 10;

	firstb >>= 18;
	secondb >>= 18;
	thirdb <<= 2;
	fourthb <<= 14;
	fifthb <<= 2;
	sixthb <<= 12;
	seventhb <<= 4;
	eighthb <<= 26;

	result |= first;
	result |= second;
	result |= third;
	result |= fourth;
	result |= fifth;
	result |= sixth;
	result |= seventh;
	result |= eighth;

	result |= firstb;
	result |= secondb;
	result |= thirdb;
	result |= fourthb;
	result |= fifthb;
	result |= sixthb;
	result |= seventhb;
	result |= eighthb;

	return result;
}

//***************************************************

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void Random::init_by_array(/*unsigned long*/uint32 init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}



//***************************************************

/* initializes mt[N] with a seed */
void Random::init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<NN; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}


//***************************************************
//********* Knuth related ***************************


//***************************************************
//  This program by D E Knuth is in the public domain and freely copyable
// Seminumerical Algorithms, 3rd edition, Section 3.6
// *    (or in the errata to the 2nd edition --- see
// *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
// *    in the changes to Volume 2 on pages 171 and following). 
// All functions below (ranf_array(), ranf_start(), ranf_arr_cycle()) are Knuth's generator.
// Tests have shown that it is extremely fast AND gives good quality data similar to the
// much slower Mersenne Twister algorithm
double Random::ran1()
{
  static double rands[NN] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 
	}; 
  static int current = NN;
  if(current == NN){
    ranf_array(rands, NN);
    current = 0;
  }
//  static int cnt = 0;
//  double temp = rands[current++];
//  cout <<cnt<<". "<<temp<<endl;
//  cnt++;
//  return(temp);
  return(rands[current++]);
}

//***************************************************

void Random::ranf_array(double aa[], int n)
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_u[j];
  for (;j<n;j++) aa[j]=mod_sum(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_u[i]=mod_sum(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_u[i]=mod_sum(aa[j-KK],ran_u[i-LL]);
}



//***************************************************

void Random::ranf_start(uint32 seed)
{
  register int t,s,j;
  double u[KK+KK-1];
  double ulp=(1.0/(1L<<30))/(1L<<22);               /* 2 to the -52 */
  double ss=2.0*ulp*((seed&0x3fffffff)+2);
//  double ss=2.0*ulp*((seed&0x000000003fffffff)+2);

//  cerr <<"Random::ranf_start seed = "<<seed<<endl;

  for (j=0;j<KK;j++) {
    u[j]=ss;                                /* bootstrap the buffer */
    ss+=ss; if (ss>=1.0) ss-=1.0-2*ulp;  /* cyclic shift of 51 bits */
  }
  u[1]+=ulp;                     /* make u[1] (and only u[1]) "odd" */
  for (s=seed&0x3fffffff,t=TT-1; t; ) {
//  for (s=seed&0x000000003fffffff,t=TT-1; t; ) {
    for (j=KK-1;j>0;j--)
      u[j+j]=u[j],u[j+j-1]=0.0;                         /* "square" */
    for (j=KK+KK-2;j>=KK;j--) {
      u[j-(KK-LL)]=mod_sum(u[j-(KK-LL)],u[j]);
      u[j-KK]=mod_sum(u[j-KK],u[j]);
    }
    if (is_odd(s)) {                             /* "multiply by z" */
      for (j=KK;j>0;j--) u[j]=u[j-1];
      u[0]=u[KK];                    /* shift the buffer cyclically */
      u[LL]=mod_sum(u[LL],u[KK]);
    }
    if (s) s>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_u[j+KK-LL]=u[j];
  for (;j<KK;j++) ran_u[j-LL]=u[j];
  for (j=0;j<10;j++) ranf_array(u,KK+KK-1);  /* warm things up */
  ranf_arr_ptr=&ranf_arr_started;
}


//***************************************************

double Random::ranf_arr_cycle()
{
  if (ranf_arr_ptr==&ranf_arr_dummy)
    ranf_start(314159L); /* the user forgot to initialize */
  ranf_array(ranf_arr_buf,QUALITY);
  ranf_arr_buf[KK]=-1;
  ranf_arr_ptr=ranf_arr_buf+1;
  return ranf_arr_buf[0];
}

//***************************************************
}
