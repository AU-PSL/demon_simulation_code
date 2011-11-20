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

class RandCache {
public:
	__m128d r;
	double l, h;
	
	RandCache(const __m128d r_ = _mm_set1_pd(0.0), 
			  const double l_ = 0.0, const double h_ = 0.0) 
	: r(r_), l(l_), h(h_) {}
};

#endif // RANDOMNUMBERS
