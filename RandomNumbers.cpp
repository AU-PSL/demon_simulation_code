//
//  RandomNumbers.cpp
//  Simulation
//
//  Created by Mark Cianciosa on 10/31/11.
//  Copyright (c) 2011 Plasma Sciences Laboratory. All rights reserved.
//

#include "RandomNumbers.h"
#include <ctime>

RandomNumbers::RandomNumbers(const double mean, const double sigma) 
: engine((unsigned int)time(NULL)), zeroToOne(0.0, 1.0), 
zeroToTwoPi(0.0, 2.0*M_PI), normalDist(mean, sigma) {}

const double RandomNumbers::uniformZeroToOne() {
	return zeroToOne(engine);
}

const double RandomNumbers::uniformZeroToTwoPi() {
	return zeroToTwoPi(engine);
}

const double RandomNumbers::guassian() {
	return normalDist(engine);
}
