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

#endif // RANDOMNUMBERS
