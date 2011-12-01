/*===- RandomNumbers.h - libSimulation -========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef RANDOMNUMBERS
#define RANDOMNUMBERS

#include <random>
#include "VectorCompatibility.h"

class RandomNumbers {
public:
	RandomNumbers();
	~RandomNumbers() {}
	
	const double uniformZeroToOne();
	const double uniformZeroToTwoPi();
	const double guassian(std::normal_distribution<double> &dist);
	
private:
	std::mt19937_64 engine;
	
	std::uniform_real_distribution<double> zeroToOne;
	std::uniform_real_distribution<double> zeroToTwoPi;
};

// Structure to hold precomuted random numbers for use with thermal forces.
struct RandCache {
	doubleV r;
    double r1, r2
#ifdef __AVX__
    r3, r4
#endif
    ;
    
	RandCache(const doubleV r_ = _mm_set1_pd(0.0), 
			  const double r1_ = 0.0, const double r2_ = 0.0
#ifdef __AVX__
              const double r3_ = 0.0, const double r4_ = 0.0 
#endif             
              ) 
	: r(r_), r1(r1_), r2(r2_) 
#ifdef __AVX__
    r3(r3_), r4(r4_) 
#endif
    {}
};

#endif // RANDOMNUMBERS
