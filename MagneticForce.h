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

class MagneticForce1D : public Force
{	
public:
	MagneticForce1D(Cloud * const myCloud);
	virtual ~MagneticForce1D() {} // destructor
	
// public functions:
	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4
	
	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);
};

class MagneticForce2D : public MagneticForce1D
{	
public:
	MagneticForce2D(Cloud * const myCloud, const double Bz);
	virtual ~MagneticForce2D() {}
	
// public functions:
	virtual void force1(const double currentTime);
	virtual void force2(const double currentTime);
	virtual void force3(const double currentTime);
	virtual void force4(const double currentTime);
	
	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

protected:
// protected variables:
	double Bz; //[T]

private:
//private functions:
	void force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentCharge);
};

class MagneticForce3D : public MagneticForce2D
{	
public:
	MagneticForce3D(Cloud * const myCloud, const double Bz, const double By, const double Bx);
	~MagneticForce3D() {}
	
// public functions:
	virtual void force1(const double currentTime);
	virtual void force2(const double currentTime);
	virtual void force3(const double currentTime);
	virtual void force4(const double currentTime);
	
	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

protected:
// protected variables:
	double By; //[T]
	double Bx; //[T]

private:
//private functions:
	void force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentVelocityZ, const __m128d currentCharge);
};

#endif // MAGNETICFORCE_H
