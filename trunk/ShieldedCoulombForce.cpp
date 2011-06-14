/*===- ShieldedCoulombForce.cpp - libSimulation -===============================
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
	for (cloud_index currentParticle = 0; currentParticle < cloud->n - 1; currentParticle += 2) 
	{
		force(currentParticle, cloud->getq1_pd(currentParticle), cloud->getEx1_pd(currentParticle), cloud->getEy1_pd(currentParticle));
	}
}

void ShieldedCoulombForce::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0; currentParticle < cloud->n - 1; currentParticle += 2) 
	{
		force(currentParticle, cloud->getq2_pd(currentParticle), cloud->getEx2_pd(currentParticle), cloud->getEy2_pd(currentParticle));
	}
}

void ShieldedCoulombForce::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0; currentParticle < cloud->n - 1; currentParticle += 2) 
	{
		force(currentParticle, cloud->getq3_pd(currentParticle), cloud->getEx3_pd(currentParticle), cloud->getEy3_pd(currentParticle));
	}
}

void ShieldedCoulombForce::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0; currentParticle < cloud->n - 1; currentParticle += 2) 
	{
		force(currentParticle, cloud->getq4_pd(currentParticle), cloud->getEx4_pd(currentParticle), cloud->getEy4_pd(currentParticle));
	}
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const __m128d charge, const __m128d Exfield, const __m128d Eyfield)
{	
	const __m128d forcevX = charge*Exfield;
	const __m128d forcevY = charge*Eyfield;

	double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);
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
