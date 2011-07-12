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
#include "VectorCompatibility.h"

class MagneticForce : public Force
{	
public:
	MagneticForce(Cloud * const myCloud, const double magneticField);
	~MagneticForce() {} // destructor
	
// public functions:
	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4
	
	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

protected:
// protected variables:
	double BField; //[T]

private:
//private functions:
	void force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentCharge);
};

#endif // MAGNETICFORCE_H
