/**
* @file  MagneticForce.cpp
* @class MagneticForce MagneticForce.h
*
* @brief Computes a magnetic force that acts in the z-direction
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "MagneticForce.h"

void MagneticForce::force1(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
		force(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle));
    END_PARALLEL_FOR
}

void MagneticForce::force2(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
		force(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle));
    END_PARALLEL_FOR
}

void MagneticForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle));
    END_PARALLEL_FOR
}

void MagneticForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
		force(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle));
    END_PARALLEL_FOR
}

/**
* @brief Computes the magnetic force with form F = q*(v x B)
*
* @param[in] currentParticle  The particle whose force is being computed
* @param[in] currentVelocityX The x-velocity of the current particle
* @param[in] currentVelocityY The y-velocity of the current particle
**/
inline void MagneticForce::force(const cloud_index currentParticle, const doubleV currentVelocityX, 
                                 const doubleV currentVelocityY) {
	const doubleV qB = mul_pd(load_pd(cloud->charge + currentParticle), BField);

	plusEqual_pd(cloud->forceX  + currentParticle, mul_pd(qB, currentVelocityY));
	minusEqual_pd(cloud->forceY + currentParticle, mul_pd(qB, currentVelocityX));
}

void MagneticForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the magnetic force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add MagneticForce bit:
		forceFlags |= MagneticForceFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), 
                            &forceFlags, const_cast<char *> ("Force configuration."), error);
	}
	
	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("magneticField"), 
                           BField, 6, const_cast<char *> ("[T] (MagneticForce)"), error);
}

void MagneticForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("magneticField"), &BField, NULL, error);
}
