/*===- RandomNumbers.cpp - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "RandomNumbers.h"
#include <chrono>

using namespace std::chrono;

RandomNumbers::RandomNumbers() 
: engine(static_cast<uint64_t> (system_clock::to_time_t(system_clock::now()))), 
zeroToOne(0.0, 1.0), zeroToTwoPi(0.0, 2.0*M_PI) {}

// Uniformly distributed random numbers between 0 - 1
const double RandomNumbers::uniformZeroToOne() {
	return zeroToOne(engine);
}

// Uniformly distributed random numbers between 0 - 2 Pi
const double RandomNumbers::uniformZeroToTwoPi() {
	return zeroToTwoPi(engine);
}

// Gaussian distributed random numbers.
const double RandomNumbers::gaussian(std::normal_distribution<double> &dist) {
	return dist(engine);
}
