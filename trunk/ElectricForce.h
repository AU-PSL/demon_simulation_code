/*===- ElectricForce.h - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
* This file is a member of the DEMON BETA project to create a more physical
* dusty plasma simulator
*
*===-----------------------------------------------------------------------===*/

#ifndef ELECTRICFORCE_H
#define ELECTRICFORCE_H

#include "Force.h"

class ElectricForce : public Force {
public:
	ElectricForce(Cloud * const C, double plasmaRad, double electricField)
	: Force(C), electric(electricField), radius(plasmaRad) {}
	~ElectricForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
	double electric, radius; // [V/m^2],[m]

	void force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY);
};

#endif // ELECTRICFORCE_H
