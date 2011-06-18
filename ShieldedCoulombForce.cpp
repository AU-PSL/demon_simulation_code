/*=== ShieldedCoulombForce.cpp - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ShieldedCoulombForce.h"
#include <cmath>

ShieldedCoulombForce::ShieldedCoulombForce(Cloud * const myCloud, const double shieldingConstant)
: Force(myCloud), shielding(shieldingConstant) {} //TODO: Shielding constant deprecated

void ShieldedCoulombForce::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n/2; currentParticle < numParticles; currentParticle++)
		force(currentParticle, cloud->getq1_pd(currentParticle*2),
			cloud->Ex[currentParticle], cloud->Ey[currentParticle]);
}

void ShieldedCoulombForce::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n/2; currentParticle < numParticles; currentParticle++)
		force(currentParticle, cloud->getq2_pd(currentParticle*2),
			  cloud->Ex[currentParticle], cloud->Ey[currentParticle]);
}

void ShieldedCoulombForce::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n/2; currentParticle < numParticles; currentParticle++)
		force(currentParticle, cloud->getq3_pd(currentParticle*2),
			  cloud->Ex[currentParticle], cloud->Ey[currentParticle]);
}

void ShieldedCoulombForce::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n/2; currentParticle < numParticles; currentParticle++)
		force(currentParticle, cloud->getq4_pd(currentParticle*2),
			  cloud->Ex[currentParticle], cloud->Ey[currentParticle]);
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const __m128d charge, const __m128d Exfield, const __m128d Eyfield)
{
	cloud->forceX[currentParticle] += charge*Exfield;
	cloud->forceY[currentParticle] += charge*Eyfield;
}


void ShieldedCoulombForce::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the drag force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ShieldedCoulombForce bit:
		forceFlags |= ShieldedCoulombForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword.
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("shieldingConstant"), shielding, 6, const_cast<char *> ("[m^-1] (ShieldedCoulombForce)"), error);
}

void ShieldedCoulombForce::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("shieldingConstant"), &shielding, NULL, error);
}
