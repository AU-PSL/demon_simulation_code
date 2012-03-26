/*===- MagneticForce.h - libSimulation -========================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef MAGNETICFORCE_H
#define MAGNETICFORCE_H

#include "Force.h"

class MagneticForce : public Force {
public:
	MagneticForce(Cloud * const C, const double magneticField)
	: Force(C), BField(magneticField) {}
	~MagneticForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4
	
	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

protected:
	double BField; // [T]

private:
	void force(const cloud_index currentParticle, const doubleV currentVelocityX, const doubleV currentVelocityY);
};

#endif // MAGNETICFORCE_H
