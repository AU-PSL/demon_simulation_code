/*===- RectConfinementForce.h - libSimulation -=================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef RECTCONFINEMENTFORCE_H
#define RECTCONFINEMENTFORCE_H

#include "Force.h"

class RectConfinementForce : public Force {
public:
	RectConfinementForce(Cloud * const C, double confineConstX, double confineConstY)
	: Force(C), confineX(confineConstX), confineY(confineConstY) {}
	// IMPORTANT: In the above constructor, confineConst_'s must be positive!
	~RectConfinementForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double confineX, confineY; // [V/m^2]

	void force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
};

#endif // RECTCONFINEMENTFORCE_H
