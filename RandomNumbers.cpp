/**
* @file  RandomNumbers.cpp
* @class RandomNumbers RandomNumbers.h
*
* @brief Provides random number generation
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
*/

#include "RandomNumbers.h"
#include <chrono>

using namespace std::chrono;

/**
* @brief Constructor for RandomNumbers class
**/
RandomNumbers::RandomNumbers() 
: engine(static_cast<uint64_t> (system_clock::to_time_t(system_clock::now()))), 
zeroToOne(0.0, 1.0), zeroToTwoPi(0.0, 2.0*M_PI) {}

/**
* @brief Uniformly distributed random numbers between 0 - 1
*
* @return Number between 0 - 1
**/
const double RandomNumbers::uniformZeroToOne() {
	return zeroToOne(engine);
}

/**
* @brief Uniformly distributed random numbers between 0 - 2*Pi
*
* @return Number between 0 - 2*pi
**/
const double RandomNumbers::uniformZeroToTwoPi() {
	return zeroToTwoPi(engine);
}


/**
* @brief Standard gaussian distributed random numbers
*
* @param[in] dist Object that produces values according to a normal distribution
* @return Number with gaussian distribution with mean 0 and std dev 1
**/
const double RandomNumbers::gaussian(std::normal_distribution<double> &dist) {
	return dist(engine);
}
