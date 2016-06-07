/**
* @file  ConfinementForce.cpp
* @class ConfinementForce ConfinementForce.h
*
* @brief Adds a force towards the center to keep the particles confined.
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "ConfinementForce.h"

void ConfinementForce::force1(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForce::force2(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

/**
* @brief Computes the confinement force with form F = -c*q*r
*
* @param[in] currentParticle  The particle whose force is being computed
* @param[in] currentPositionX The x-position of the current particle
* @param[in] currentPositionY The y-position of the current particle
**/
inline void ConfinementForce::force(const cloud_index currentParticle, const doubleV currentPositionX, 
                                    const doubleV currentPositionY) {
	const doubleV cV = mul_pd(load_pd(cloud->charge + currentParticle), confine);
	
	plusEqual_pd(cloud->forceX + currentParticle, mul_pd(cV, currentPositionX));
	plusEqual_pd(cloud->forceY + currentParticle, mul_pd(cV, currentPositionY));
}

void ConfinementForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the confinement force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ConfinementForce bit:
		forceFlags |= ConfinementForceFlag;
		
		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConst"), confine, 
                           6, const_cast<char *> ("[V/m^2] (ConfinementForce)"), error);
}

void ConfinementForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConst"), &confine, NULL, error);
}
