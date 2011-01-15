/*===- DragForce.h - libSimulation -============================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef DRAGFORCE_H
#define DRAGFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class DragForce1D : public Force
{
public:
	DragForce1D(Cloud * const myCloud, const double gamma);
	~DragForce1D() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); //rk substep 1
	virtual void force2(const double currentTime); //rk substep 2
	virtual void force3(const double currentTime); //rk substep 3
	virtual void force4(const double currentTime); //rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

protected:
//protected variables:
	double dragConst; //[s^-1]

private:
//private functions:
	void force(const unsigned int currentParticle, const __m128d currentVelocityX);
};

class DragForce2D : public DragForce1D
{
public:
	DragForce2D(Cloud * const myCloud, const double gamma);
	~DragForce2D() {}

private:
//private functions:
	void force(const unsigned int currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY);
};

class DragForce3D : public DragForce2D
{
public:
	DragForce3D(Cloud * const myCloud, const double gamma);
	~DragForce3D() {}

private:
//private functions:
	void force(const unsigned int currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentVelocityZ);
};

#endif /* DRAGFORCE_H */
