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
#include "VectorCompatibility.h"

class RectConfinementForce1D : public Force
{
public:
	RectConfinementForce1D(Cloud * const myCloud, double confineConstX); //confinement consts must be positive!
	~RectConfinementForce1D() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//private variables:
	double confineX;

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX);
};

class RectConfinementForce2D : public RectConfinementForce1D
{
public:
	RectConfinementForce2D(Cloud * const myCloud, double confineConstX, double confineConstY);
	~RectConfinementForce2D() {}

private:
//private variables:
	double confineY;

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
};

class RectConfinementForce3D : public RectConfinementForce2D
{
public:
	RectConfinementForce3D(Cloud * const myCloud, double confineConstX, double confineConstY, double confineConstZ);
	~RectConfinementForce3D() {}

private:
//private variables:
	double confineZ;

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d currentPositionZ);
};

#endif /* RECTCONFINEMENTFORCE_H */
