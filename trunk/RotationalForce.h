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

class RotationalForce : public Force {
public:
	RotationalForce(Cloud * const C, const double rmin, const double rmax, const double rotConst)
	: Force(C), innerRad(rmin), outerRad(rmax), rotationalConst(rotConst) {}
	~RotationalForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

private:
	double innerRad, outerRad, rotationalConst; // [m], [m], [N]

	void force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY);
    const doubleV rotationConstant(const int mask);
};

#endif // ROTATIONALFORCE_H
