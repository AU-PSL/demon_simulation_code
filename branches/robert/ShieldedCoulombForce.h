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

class ShieldedCoulombForce : public Force
{	
public:
	ShieldedCoulombForce(Cloud * const myCloud, const double shieldingConstant);	//overloaded constructor
	~ShieldedCoulombForce() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1_1D(const double currentTime); //rk substep 1
	void force2_1D(const double currentTime); //rk substep 2
	void force3_1D(const double currentTime); //rk substep 3
	void force4_1D(const double currentTime); //rk substep 4

	void force1_2D(const double currentTime); 
	void force2_2D(const double currentTime); 
	void force3_2D(const double currentTime); 
	void force4_2D(const double currentTime); 

	void force1_3D(const double currentTime); 
	void force2_3D(const double currentTime); 
	void force3_3D(const double currentTime); 
	void force4_3D(const double currentTime); 

	void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//public variables:
	double shielding;

//private functions:
	void force1D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX);
	void force2D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX, const double displacementY);
	void force3D(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX, const double displacementY, const double displacementZ);

	void force1D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX);
	void force2D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY);
	void force3D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);

	void forcer1D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX);
	void forcer2D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY);
	void forcer3D(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
};

#endif /* SHIELDEDCOULOMBFORCE_H */
