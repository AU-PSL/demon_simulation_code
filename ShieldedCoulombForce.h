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

class ShieldedCoulombForce : public Force {
public:
	ShieldedCoulombForce(Cloud * const C, const double shieldingConstant);
	~ShieldedCoulombForce();

// public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
// public variables:
	double shielding; // [m^-1]
    SEMAPHORES
	
	static const double coulomb;

// private functions:
	void force(const cloud_index currentParticle, const cloud_index iParticle, 
	           const double currentCharge, const double iCharge,
	           const double displacementX, const double displacementY);
	void force(const cloud_index currentParticle, const cloud_index iParticle, 
	           const __m128d currentCharge, const __m128d iCharge,
	           const __m128d displacementX, const __m128d displacementY);
	void forcer(const cloud_index currentParticle, const cloud_index iParticle,
	            const __m128d currentCharge, const __m128d iCharge,
	            const __m128d displacementX, const __m128d displacementY);
};

#endif // SHIELDEDCOULOMBFORCE_H
