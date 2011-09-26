//////////////////////////////////////////////////////////////////////////
//                                                                      //
// This file is distributed as part of the genomeSIM source code package//
// and may not be redistributed in any form without written permission  //
// from Dr. Marylyn Ritchie (ritchie@chgr.mc.vanderbilt.edu).           //
// Permission is granted to modify this file for your own personal      //
// use, but modified versions must retain this notice and must not be   //
// distributed.                                                         //
//                                                                      //
// This application is provided "as is" without express or implied      //
// warranty.                                                            //
//                                                                      //  
//////////////////////////////////////////////////////////////////////////

#include "random.h"


namespace Utility {


Random Random::globalGenerator=Random(1371);


double Random::Seed(long randomSeed){
	RND_LOCK
		seed = randomSeed;
		gen.RandomInit((int)randomSeed);
	RND_UNLOCK;

	return 0.0;
}

long Random::Reset(long inc) {
	RND_LOCK
		seed+=inc;
		gen.RandomInit((int)seed);
	RND_UNLOCK
	return seed;
}
/**
 * Return a floating point random number
 */
double Random::drand() { 	
	RND_LOCK
		double r = gen.Random();
	RND_UNLOCK;
//	cout<<"drand()="<<r<<"\n";

	return r;
}
/**
 * Return a long random number
 */
long Random::lrand(int min, int max) {
	if (min == max)
		return min;
	RND_LOCK;
		int r = gen.IRandom(min, max);
	RND_UNLOCK;

//	cout<<"lrand("<<min<<","<<max<<")="<<r<<"\n";
	return (long)r;
}

long Random::lrand() { 		
	return lrand(0, RAND_MAX);
}

float Random::operator()() {
	RND_LOCK;
	float v=1.0;
	while (v==1.0)
		v = drand();
	
	RND_UNLOCK;
	return v;
}
float Random::operator()(float n) {
	RND_LOCK;
	float v = 0;
	if (n > 0)
		v = drand()*n;
	RND_UNLOCK;


	return v;
}

int Random::operator()(int n) {
	int v=0;
	if (n>1)
		v = (int)lrand(0, (long)n - 1);
	return v;
}

double Random::operator()(double n) {
	RND_LOCK;
	double v = 0;
	if (n > 0.0)
		v = drand()*n;
	RND_UNLOCK;

	return v;
}

long Random::operator()(long n) {
	long v = 0;
	if (n > 0)
		v = lrand(0,n - 1);
	return v;
}

}
