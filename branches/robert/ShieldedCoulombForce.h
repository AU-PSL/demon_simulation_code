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
	virtual ~ShieldedCoulombForce1D() {}; //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); //rk substep 1
	virtual void force2(const double currentTime); //rk substep 2
	virtual void force3(const double currentTime); //rk substep 3
	virtual void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

protected:
//protected variables:
	double shielding;

private:
//private functions:
	void force(const cloud_index currentParticle, const cloud_index iParticle, const double displacementX);
	void force(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX);
	void forcer(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX);
};

class ShieldedCoulombForce2D : public ShieldedCoulombForce1D
{
public:
	ShieldedCoulombForce2D(Cloud * const myCloud, const double shieldingConstant);
	virtual ~ShieldedCoulombForce2D() {};

//public functions:
	virtual void force1(const double currentTime);
	virtual void force2(const double currentTime);
	virtual void force3(const double currentTime);
	virtual void force4(const double currentTime);

private:
//private functions:
	void force(const cloud_index currentParticle, const cloud_index iParticle, const double displacementX, const double displacementY);
	void force(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX, const __m128d displacementY);
	void forcer(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX, const __m128d displacementY);
};

class ShieldedCoulombForce3D : public ShieldedCoulombForce2D
{
public:
	ShieldedCoulombForce3D(Cloud * const myCloud, const double shieldingConstant);
	~ShieldedCoulombForce3D() {};

//public functions:
	void force1(const double currentTime);
	void force2(const double currentTime);
	void force3(const double currentTime);
	void force4(const double currentTime);

private:
//private functions:
	void force(const cloud_index currentParticle, const cloud_index iParticle, const double displacementX, const double displacementY, const double displacementZ);
	void force(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
	void forcer(const cloud_index currentParticle, const cloud_index iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
};

#endif /* SHIELDEDCOULOMBFORCE_H */
