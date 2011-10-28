/*===- RectConfinementForce.cpp - libSimulation -===============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

//TODO: Charge-dependence

#include "RectConfinementForce.h"
#include "Parallel.h"
	
RectConfinementForce::RectConfinementForce(Cloud * const myCloud, double confineConstX, double confineConstY)
: Force(myCloud), confineX(-confineConstX), confineY(-confineConstY) {}

void RectConfinementForce::force1(const double currentTime) {
    (void)currentTime;
	begin_parallel_for(currentParticle, numParticles, cloud->n, 2)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    end_parallel_for
}

void RectConfinementForce::force2(const double currentTime) {
    (void)currentTime;
	begin_parallel_for(currentParticle, numParticles, cloud->n, 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    end_parallel_for
}

void RectConfinementForce::force3(const double currentTime) {
    (void)currentTime;
	begin_parallel_for(currentParticle, numParticles, cloud->n, 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    end_parallel_for
}

void RectConfinementForce::force4(const double currentTime) {
    (void)currentTime;
	begin_parallel_for(currentParticle, numParticles, cloud->n, 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    end_parallel_for
}

inline void RectConfinementForce::force(const cloud_index currentParticle, const __m128d currentPositionX, 
                                        const __m128d currentPositionY) {
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_set1_pd(confineX)*currentPositionX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + _mm_set1_pd(confineY)*currentPositionY);
}

void RectConfinementForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the rectangular confinement force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add RectConfinementForce bit:
		forceFlags |= RectConfinementForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConstX"), confineX, 
                           6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("confineConstY"), confineY, 
                           6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
	}
}

void RectConfinementForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error) {
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConstX"), &confineX, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("confineConstY"), &confineY, NULL, error);
	}
}
