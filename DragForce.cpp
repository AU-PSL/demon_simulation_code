/*===- DragForce.cpp - libSimulation -==========================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "DragForce.h"

void DragForce::force1(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle));
    END_PARALLEL_FOR
}

void DragForce::force2(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle));
    END_PARALLEL_FOR
}

void DragForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle));
    END_PARALLEL_FOR
}

void DragForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle));
    END_PARALLEL_FOR
}

inline void DragForce::force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY) {
	const __m128d drag = _mm_set1_pd(dragConst)*_mm_load_pd(cloud->mass + currentParticle);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + drag*currentVelocityX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + drag*currentVelocityY);
}

void DragForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the drag force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add DragForce bit:
		forceFlags |= DragForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}
	
	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("dragConst"), dragConst, 
                           6, const_cast<char *> ("[s^-1] (DragForce)"), error);
}

void DragForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("dragConst"), &dragConst, NULL, error);
}
