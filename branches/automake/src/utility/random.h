// Ran2.h

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

///////////////////////////////////////////////////////////////////// 
//
// Pseudo-random number generator used within the library.
// Uses standard C/C++ random number generator.
// Can be replaced easily by changing these 
//
/////////////////////////////////////////////////////////////////////

#ifndef __RANDOM_NUMBER_GENERATOR_H__
#define __RANDOM_NUMBER_GENERATOR_H__
#include <iostream>
#include "types.h"
//using namespace std;


#ifdef USE_MPI
#include <pthread.h>
#define RND_INIT pthread_mutex_init(&lock, NULL);
#define RND_DESTROY pthread_mutex_destroy(&lock);

#define RND_LOCK pthread_mutex_lock(&lock);
#define RND_UNLOCK pthread_mutex_unlock(&lock);

#else
//We don't want these to mean anything if we arne't using MPI

#define RND_INIT
#define RND_DESTROY
#define RND_LOCK
#define RND_UNLOCK

#endif 	//USE_MPI

#include "random/randomc.h"

namespace Utility {

using namespace std;

typedef CRandomMersenne TRandomMersenne;

/**
 * @brief Base class for all random number generators. 
 * @Note It would seem reasonable to allow the application (or configuration) to set this object up and 
 * pass it to any object that requires a random number generator 
 */
class Random {

public:
	Random(long seed=1371) : gen(seed) {
		RND_INIT
	}
	virtual ~Random() {
		RND_DESTROY
	}

	virtual double Seed(long randomSeed);
	virtual long lrand(); 
	virtual long lrand(int min, int max);

	virtual double drand();

	virtual long Reset(long increment = 0);

	/**
	 * Used by STL
	 */
	long operator()(long v);
	int operator()(int v);
	double operator()(double n);
	float operator()(float v);
	float operator()();



	/**
	 * Keep up with the seed used
	 */
	long GetSeed() { 			
		return seed; 
	}

	static Random globalGenerator;

protected:
#ifdef USE_MPI
	///<Since the calculation is mildly complex and we want deterministic results for a given seed, we should be safe
	pthread_mutex_t lock;			
#endif
	TRandomMersenne gen;
	long seed;
};




}

#endif

