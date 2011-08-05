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
#include "VectorCompatibility.h"

//ConfinementForceVoid1D inherits from ConfinementForce3D (as opposed to 
// ConfinementForce1D) to circumvent the need for multiple inheritence
// in ConfinementForce2D/3D.
class ConfinementForceVoid1D : public ConfinementForce3D
{	
public:
	ConfinementForceVoid1D(Cloud * const myCloud, double confineConst, double plasmaPotential, double voidDecay);
	// IMPORTANT: In the above constructor, confineConst must be positive!
	virtual ~ConfinementForceVoid1D() {} // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

protected:
// protected variables:
	double confine, decay, potentialOffset;

private:
// private functions:
	void force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d charge);	// common force calculator
};

class ConfinementForceVoid2D : public ConfinementForceVoid1D
{	
public:
	ConfinementForceVoid2D(Cloud * const myCloud, double confineConst, double plasmaPotential, double voidDecay);
	// IMPORTANT: In the above constructor, confineConst must be positive!
	virtual ~ConfinementForceVoid2D() {} // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

private:
// private functions:
	void force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d charge);	// common force calculator
};

class ConfinementForceVoid3D : public ConfinementForceVoid2D
{	
public:
	ConfinementForceVoid3D(Cloud * const myCloud, double confineConst, double plasmaPotential, double voidDecay);
	// IMPORTANT: In the above constructor, confineConst must be positive!
	virtual ~ConfinementForceVoid3D() {} // destructor
};

#endif // CONFINEMENTFORCEVOID_H
