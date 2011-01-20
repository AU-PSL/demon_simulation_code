/*===- RotationalForce.h - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef ROTATIONALFORCE_H
#define ROTATIONALFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class RotationalForce : public Force
{	
public:
	RotationalForce(Cloud * const myCloud, const double rmin, const double rmax, const double rotConst);
	~RotationalForce() {} // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

private:
// private variables:
	double innerRad;
	double outerRad;
	double rotationalConst;

// private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
};

#endif // ROTATIONALFORCE_H
