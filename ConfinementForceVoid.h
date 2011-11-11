/*===- ConfinementForceVoid.h - libSimulation -=================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef CONFINEMENTFORCEVOID_H
#define CONFINEMENTFORCEVOID_H

#include "ConfinementForce.h"

class ConfinementForceVoid : public ConfinementForce {
public:
	ConfinementForceVoid(Cloud * const C, double confineConst, double voidDecay)
	: ConfinementForce(C, confineConst), decay(voidDecay) {}
	
	// IMPORTANT: In the above constructor, confineConst must be positive!
	~ConfinementForceVoid() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double confine, decay; // [V/m^2], [m^-1]

	void force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);	// common force calculator
};

#endif // CONFINEMENTFORCEVOID_H
