/*===- ShieldedCoulombForce.h - libSimulation -=================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef SHIELDEDCOULOMBFORCE_H
#define SHIELDEDCOULOMBFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class ShieldedCoulombForce1D : public Force
{
public:
	ShieldedCoulombForce1D(Cloud * const myCloud, const double shieldingConstant);
	~ShieldedCoulombForce1D() {}; //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
//public variables:
	double shielding;

//private functions:
	void force1D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX);
	void force1D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX);
	void forcer1D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX);
};

class ShieldedCoulombForce2D : public ShieldedCoulombForce1D
{
public:
	ShieldedCoulombForce2D(Cloud * const myCloud, const double shieldingConstant);
	~ShieldedCoulombForce2D() {};

private:
//private functions:
	void force2D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX, const double displacementY);
	void force2D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY);
	void forcer2D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY);
};

class ShieldedCoulombForce3D : public ShieldedCoulombForce2D
{
public:
	ShieldedCoulombForce3D(Cloud * const myCloud, const double shieldingConstant);
	~ShieldedCoulombForce2D() {};

private:
//private functions:
	void force3D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX, const double displacementY, const double displacementZ);
	void force3D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
	void forcer3D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
};

#endif /* SHIELDEDCOULOMBFORCE_H */
