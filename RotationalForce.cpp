/**
* @file  RotationalForce.cpp
* @class RotationalForce RotationalForce.h
*
* @brief Computes a rotational force in the azimuthal direction
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "RotationalForce.h"

void RotationalForce::force1(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force2(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void RotationalForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}


/**
* @brief Computes a rotational force if a particle is within a radial range
*        with form F_x = -c*y/r. F_y = c*x/r
*
* @param[in] currentParticle  The particle whose force is being computed
* @param[in] currentPositionX The x-position of the current particle
* @param[in] currentPositionY The y-position of the current particle
**/
inline void RotationalForce::force(const cloud_index currentParticle, 
                                   const doubleV currentPositionX, 
                                   const doubleV currentPositionY) {
    const doubleV dustRadV = length_pd(currentPositionX, currentPositionY);

	// dustRad > innerRad && dustRadV < outerRad
	const int mask = movemask_pd(and_pd(cmpgt_pd(dustRadV, innerRad), 
                                        cmplt_pd(dustRadV, outerRad)));
	if (!mask)
		return; // niether in, early return
	
	doubleV cRotConst = select_pd(mask, rotationalConst, 0.0);
	
	// force in theta direction:
	minusEqual_pd(cloud->forceX + currentParticle, div_pd(mul_pd(cRotConst, currentPositionY), dustRadV));
	plusEqual_pd(cloud->forceY + currentParticle, div_pd(mul_pd(cRotConst, currentPositionX), dustRadV));
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
