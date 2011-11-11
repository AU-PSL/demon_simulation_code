/*===- RotationalForce.cpp - libSimulation -====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "RotationalForce.h"

void RotationalForce::force1(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force2(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

inline void RotationalForce::force(const cloud_index currentParticle, const __m128d currentPositionX, 
                                   const __m128d currentPositionY) {
	const __m128d dustRadV = _mm_sqrt_pd(currentPositionX*currentPositionX + currentPositionY*currentPositionY);

	// dustRad > innerRad && dustRadV < outerRad
	const int mask = _mm_movemask_pd(_mm_and_pd(_mm_cmpgt_pd(dustRadV, _mm_set1_pd(innerRad)), 
												_mm_cmplt_pd(dustRadV, _mm_set1_pd(outerRad))));
	if (!mask)
		return; // niether in, early return
	
	__m128d cRotConst = _mm_set_pd((mask & 2) ? rotationalConst : 0.0, // _mm_set_pd() is backwards.
								   (mask & 1) ? rotationalConst : 0.0);
	
	// force in theta direction:
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	
	// Fx = -c*x/r;
	// Fy = c*y/r;
	_mm_store_pd(pFx, _mm_load_pd(pFx) - cRotConst*currentPositionY/dustRadV);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + cRotConst*currentPositionX/dustRadV);
}

void RotationalForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the drag force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add RotationalForce bit:
		forceFlags |= RotationalForceFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("rotationalConst"), rotationalConst, 
                           6, const_cast<char *> ("[N] (RotationalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("innerRadius"), innerRad, 
                           6, const_cast<char *> ("[m] (RotationalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("outerRadius"), outerRad, 
                           6, const_cast<char *> ("[m] (RotationalForce)"), error);
	}
}

void RotationalForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error) {
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("rotationalConst"), &rotationalConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("innerRadius"), &innerRad, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("outerRadius"), &outerRad, NULL, error);
	}
}
