#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

namespace Simla {

using namespace std;


//typedef   signed long int32;     
//typedef unsigned long uint32;     
typedef   signed int int32;     
typedef unsigned int uint32;    

#if defined (WIN32)
	#include <sys/timeb.h>
	#include <time.h>
	#include <process.h>
    #include <string.h>
	#include <io.h>
//	#include <strstrea.h>
#endif

#if defined (SOLARIS)  || defined (OSX)
	#include <sys/time.h>
	#include <sys/timeb.h>
	#include <unistd.h>
	#include <sstream>
	#include <backward/strstream>
	using std::ostrstream;
#endif



#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
/* Period parameters */  
#define NN 624

#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

// related to Knuth's new random number generator 
#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define mod_sum(x,y) (((x)+(y))-(int)((x)+(y)))   /* (x+y) mod 1.0 */
#define QUALITY 1009 /* recommended quality level for high-res use */
#define TT  70   /* guaranteed separation between streams */
#define is_odd(s) ((s)&1)
#define ranf_arr_next() (*ranf_arr_ptr>=0? *ranf_arr_ptr++: ranf_arr_cycle())

//  This program by D E Knuth is in the public domain and freely copyable
// Seminumerical Algorithms, 3rd edition, Section 3.6
// *    (or in the errata to the 2nd edition --- see
// *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
// *    in the changes to Volume 2 on pages 171 and following). 
// end Knuth

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

const long int a1 = 16807;
const long int b15 = 32768L;
const long int b16 = 65536L;
const long int p = 2147483647L;
//const uint32 maxint = 2147483647;
const double sqrt2pi = 2.50663;
const int def_init = 4;

class Random{

 
 public:
  Random();
  ~Random();
  // It's really important that we only have one instance of these functions.
  // We need to be able to produce the same sequence of random numbers even if
  // we have several instances of the class. Like one in a.cpp, b.cpp & c.cpp
  // They all need to utilize the SAME randomness function to secure reproducible results.
  static void setseed(uint32 i);                         // inits "idum" and sets static variable seed to 1
  static void clearseed();                                     // sets static variable seed to 0
  static void init();
  static double ran1();                                // same algorithm, just return int
  static uint32 ran2();
//  static double rana();

  static double ranb();
  static double ranc();
  static void ranf_array(double aa[], int n);
  static void ranf_start(uint32 seed);
  static double ranf_arr_cycle();

  static double rnorm(double mu = 0, double dev = 1);
  static double normdist(double mu, double dev, double x);      // given mu/(=mean), standard deviation and x, returns f(x) where f is normal distribution with mu & std.dev.
  static int ranint(int);                                      // given arg = n > 0, returns an integer from 0,1,...,(n-1) with uniform distribution
  static void sgenrand(unsigned long);
  static int maxbits;			// bits in type "unsigned long"
  static uint32 maxint;			// largest unsigned int
  static int negmaxbits;		// (-1) * maxbits
  static uint32 shuffle_seed(uint32);
  static uint32 mt[NN];
  uint32 testseed();                                     // to test variability of seed generation
  static void bitsinunsignedlong();								// count number of bits in uint32
  static void init_by_array(/*unsigned long*/uint32 init_key[], int key_length = def_init);
  static void init_genrand(unsigned long);
  static void bye();                                            // writes new ".randosa2" file with 4 seeds for the next run.
 private:
	 // The all important seed value. 
	 // seed == 0 => be as random as possible, we don't want to reproduce the result (see rana())
	 // seed != 0 => we want the exact same sequence of random numbers every time we
	 // start the generator with that particular seed value (see ranb())


};
extern Random r;

}
#endif
