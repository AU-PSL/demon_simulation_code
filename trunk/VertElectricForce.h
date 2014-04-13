/*===- VertElectricForce.h - libSimulation -=======================================
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

#ifndef VERTELECTRICFORCE_H
#define VERTELECTRICFORCE_H

#include "Force.h"

class VertElectricForce : public Force {
public:
	VertElectricForce(Cloud * const C, double vertElectricField, double vertDecay)
	: Force(C), vertElectric(vertElectricField), vertDec(vertDecay) {}
	~VertElectricForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
	double vertElectric, vertDec; // [V/m^2],[m],[m]

	void force(const cloud_index currentParticle, const doubleV currentPositionY);
};

#endif // VERTELECTRICFORCE_H
