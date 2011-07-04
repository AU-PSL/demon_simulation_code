/*===- ConfinementForce.cpp - libSimulation -===================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ConfinementForce.h"
	
ConfinementForce::ConfinementForce(Cloud * const myCloud, double confineConst, double plasmaPotential) : Force(myCloud), confine(confineConst), potentialOffset(plasmaPotential) {}

void ConfinementForce::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void ConfinementForce::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void ConfinementForce::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void ConfinementForce::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

inline void ConfinementForce::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d charge)
{
	const __m128d cV = _mm_set1_pd(confine);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + cV*charge*currentPositionX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + cV*charge*currentPositionY);

	const __m128d rr = currentPositionX*currentPositionX + currentPositionY*currentPositionY;
	double * pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_set1_pd(-0.5)*cV*rr + _mm_set1_pd(potentialOffset));
}

void ConfinementForce::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the confinement force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ConfinementForce bit:
		forceFlags |= ConfinementForceFlag;	// compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConst"), confine, 6, const_cast<char *> ("[V/m^2] (ConfinementForce)"), error);

	// write background plasma potential offset:
	if(!*error)
		fits_write_key_dbl(file, const_cast<char *> ("plasmaPotential"), potentialOffset, 6, const_cast<char *> ("[V] (background plasma potential offset)"), error);
}

void ConfinementForce::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConst"), &confine, NULL, error);
}
