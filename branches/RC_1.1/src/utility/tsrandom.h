#ifndef UTILITY_TSRANDOM_H
#define UTILITY_TSRANDOM_H

/**
 * Creates a thread-safe singleton random generator engine compatible with both
 * the Boost and STL generator requirements.
 *
 * This class wraps a given generator (by default the boost Mersenne Twister)
 * into a thread-safe singleton class.  The generator that is given MUST comply
 * with the boost.Random 1.47 requirements for Uniform Number Generators.
 *
 * Note that since we are assuming the use of the Boost.Random library, we will
 * also assume the use of the Boost.Threading library for the mutex locking
 */

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#ifdef USE_MPI
#include <boost/thread/mutex.hpp>
#define RND_LOCK mx.lock();
#define RND_UNLOCK mx.unlock();
#else
//We don't want these to mean anything if we arne't using MPI

#define RND_LOCK
#define RND_UNLOCK

#endif 	//USE_MPI

namespace Utility{

template<class T = boost::random::mt19937>
class TSRandom {

public:
	typedef typename T::result_type result_type;

	TSRandom() {
		InitGenerator();
	}
	explicit TSRandom(int S) {
		InitGenerator(S);
	}

	/**
	 * Operator()() that allows for usage by Boost.
	 * This function is equivalent to the operator()() given by the PRNG used.
	 *
	 * @return A T::return_type distributed according to the given generator.
	 */
	result_type operator()() {
		result_type val;
		RND_LOCK
		val = (*gen)();
		RND_UNLOCK
		return val;
	}

	/**
	 * Operator()(long) that allows for usage by STL
	 *
	 * @return A long uniformly distributed in the range [0, N)
	 */
	long operator()(long N){
		long val;
		RND_LOCK
		val = unif01_dist(*gen)*N;
		RND_UNLOCK
		return val;
	}

	/**
	 * Sets the seed of the PRNG.
	 *
	 * @param An integer seed for the PRNG
	 */
	void seed(int S){
		RND_LOCK
		gen->seed(S);
		RND_UNLOCK
	}

	/**
	 * Maximum value of the generator
	 *
	 * @return The maximum value of the generator used.
	 */
	result_type max(){ return gen->max(); }
	/**
	 * Minimum value of the generator
	 *
	 * @return The minimum value of the generator used.
	 */
	result_type min(){ return gen->min(); }

	/**
	 * Returns a value uniformly distributed in [0,1).
	 *
	 * Provided for backwards compatibility
	 *
	 * @return A sample from Unif([0,1))
	 */
	double drand(){
		double val;
		RND_LOCK
		val = unif01_dist(*gen);
		RND_UNLOCK
		return val;
	}

	/**
	 * Returns a reference to a predefined global generator.
	 */
	static TSRandom<T>& globalGenerator(){
		return *(global_gen ? global_gen : global_gen = new TSRandom<T>());
	}

private:
	void InitGenerator() {
		RND_LOCK
		if (gen == 0) {
			gen = new T();
		}
		RND_UNLOCK
	}

	void InitGenerator(int S) {
		RND_LOCK
		if (gen == 0) {
			gen = new T(S);
		}
		RND_UNLOCK
	}

	//! A Helper distribution to generate #s in [0,1)
	static boost::random::uniform_01<> unif01_dist;
	//! The actual PRNG
	static T* gen;
	//! A pointer to an instance of this class (for speed / backward compatibility)
	static TSRandom<T>* global_gen;
#ifdef USE_MPI
	///<Since the calculation is mildly complex and we want deterministic results for a given seed, we should be safe
	static boost::mutex mx;
#endif
};

#ifdef USE_MPI
template<class T>
boost::mutex TSRandom<T>::mx;
#endif

template<class T>
boost::random::uniform_01<> TSRandom<T>::unif01_dist;

template<class T>
T* TSRandom<T>::gen = 0;

template<class T>
TSRandom<T>* TSRandom<T>::global_gen = 0;

}; //namespace Utility

#endif
