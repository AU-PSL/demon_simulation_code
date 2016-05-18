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
	const double gaussian(std::normal_distribution<double> &dist);
	
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
    
    RandCache(RandomNumbers &rands) :
#ifdef __AVX__
    r(_mm256_set_pd(rands.uniformZeroToOne(),
                    rands.uniformZeroToOne(),
                    rands.uniformZeroToOne(),
                    rands.uniformZeroToOne())),
    r1(rands.uniformZeroToTwoPi()), r2(rands.uniformZeroToTwoPi()),
    r3(rands.uniformZeroToTwoPi()), r4(rands.uniformZeroToTwoPi()) {}
#else
    r(_mm_set_pd(rands.uniformZeroToOne(),
                 rands.uniformZeroToOne())),
    r1(rands.uniformZeroToTwoPi()), r2(rands.uniformZeroToTwoPi()) {}
#endif
	RandCache() 
	: r(set1_pd(0.0)), r1(0.0), r2(0.0) 
#ifdef __AVX__
    r3(0.0), r4(0.0) 
#endif
    {}
};

#endif // RANDOMNUMBERS
