/*===- GravitationalForce.cpp - libSimulation -=====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
* This file is a member of the DEMON BETA project to create a more physical
* dusty plasma simulator
*
*===-----------------------------------------------------------------------===*/

#include "GravitationalForce.h"

void GravitationalForce::force1(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle);
    END_PARALLEL_FOR
}

void GravitationalForce::force2(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle);
    END_PARALLEL_FOR
}

void GravitationalForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle);
    END_PARALLEL_FOR
}

void GravitationalForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle);
    END_PARALLEL_FOR
}

// F = -m*g
inline void GravitationalForce::force(const cloud_index currentParticle) {
	const doubleV gravity = mul_pd(load_pd(cloud->mass + currentParticle), gravitational);

	plusEqual_pd(cloud->forceY + currentParticle, gravity);
}

void GravitationalForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the gravitational force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add GravitationalForce bit:
		forceFlags |= GravitationalForceFlag;
		
		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("GravitationalFieldStrength"), gravitational, 
                           6, const_cast<char *> ("[m/s^2] (GravitationalForce)"), error);

}

void GravitationalForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("gravitationalFieldStrength"), &gravitational, NULL, error);

}
