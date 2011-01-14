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
#include "VectorCompatibility.h"

class ConfinementForce1D : public Force
{
public:
	ConfinementForce1D(Cloud * const myCloud, double confineConst);
	//IMPORTANT: In the above constructor, confineConst must be positive!
	~ConfinementForce1D() {}

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
//private variables:
	double confine;

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX); //common force calculator
};

class ConfinementForce2D : public ConfinementForce1D
{
public:
	ConfinementForce2D(Cloud * const myCloud, double confineConst); //overloaded constructor
	~ConfinementForce2D() {} //destructor

private:
//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
};

class ConfinementForce3D : public ConfinementForce2D
{
public:
	ConfinementForce3D(Cloud * const myCloud, double confineConst); //overloaded constructor
	~ConfinementForce3D() {} //destructor

private:
//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d currentPositionZ);
};

#endif /* CONFINEMENTFORCE_H */
