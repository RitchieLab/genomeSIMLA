#ifndef __UTILITY_RANDOM_H
#define __UTILITY_RANDOM_H

#include <boost/random.hpp>

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

#endif   //USE_MPI

namespace Utility{

using namespace boost::random;

class Random: public mt19937{
  
public:
  Random() : mt19937() {RND_INIT}
  explicit Random (uint32_t seed) : mt19937(seed) {RND_INIT}
  
  ~Random() {RND_DESTROY}
    
  inline void Seed(long randomSeed){
		RND_LOCK
		this->seed(randomSeed);
		seed_val = randomSeed;
		RND_UNLOCK
  }
  
  inline long lrand(){ return lrand(this->min(), this->max()); }
  
  inline long lrand(int min_v, int max_v){ 
    uniform_int_distribution<> dist;
    RND_LOCK
    long result = dist(*this);
    RND_UNLOCK
    return result;
  }

  inline double drand(){
    uniform_01<> dist;
  	RND_LOCK
  	double result = dist(*this);
  	RND_UNLOCK
  	return result;
  }

  inline long Reset(long increment = 0){
  	RND_LOCK
		seed_val+=increment;
		this->seed(seed_val);
		RND_UNLOCK
	  return seed_val;
	}
  
  	/**
	 * Used by STL
	 */
	inline long operator()(long v){ return v>0 ? lrand(0, v-1) : 0;};
	inline int operator()(int v){return v>0 ? lrand(0, v-1) : 0;};
	inline double operator()(double v){return v>0 ? drand()*v : 0;};
	inline float operator()(float v){return v>0 ? drand()*v : 0;};
	//long operator()(){return drand();};
	inline uint32_t operator()(){return lrand();}
	
  inline long GetSeed() { 			
		return seed_val; 
	}
    
  static Random globalGenerator;
  
protected:
#ifdef USE_MPI
  ///<Since the calculation is mildly complex and we want deterministic results for a given seed, we should be safe
  pthread_mutex_t lock;      
#endif
	long seed_val;
};

}


#endif
