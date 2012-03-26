/*===- ConfinementForce.h - libSimulation -=====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef CONFINEMENTFORCE_H
#define CONFINEMENTFORCE_H

#include "Force.h"

class ConfinementForce : public Force {
public:
	ConfinementForce(Cloud * const C, double confineConst)
	: Force(C), confine(confineConst) {}
	// IMPORTANT: In the above constructor, confineConst must be positive!
	~ConfinementForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
	double confine; // [V/m^2]

	void force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY);
};

#endif // CONFINEMENTFORCE_H
