/*===- ConfinementForceVoid.cpp - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ConfinementForceVoid.h"
#include <cmath>
	
ConfinementForceVoid::ConfinementForceVoid(Cloud * const myCloud, double confineConst, double voidDecay) : 
ConfinementForce(myCloud, confineConst), decay(voidDecay) {}

void ConfinementForceVoid::force1(const double currentTime) {
	ConfinementForce::force1(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force2(const double currentTime) {
	ConfinementForce::force2(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force3(const double currentTime) {
	ConfinementForce::force3(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force4(const double currentTime) {
	ConfinementForce::force4(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

inline void ConfinementForceVoid::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY) {
	const __m128d decayV = _mm_set1_pd(decay)*_mm_load_pd(cloud->charge + currentParticle);
	const __m128d r = _mm_sqrt_pd(currentPositionX*currentPositionX + currentPositionY*currentPositionY);

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	double rL, rH;
	_mm_storel_pd(&rL, r);
	_mm_storeh_pd(&rH, r);
	const __m128d expR = _mm_set_pd(exp(-decay*rH), exp(-decay*rL));

	double xL, xH, yL, yH;
	_mm_storel_pd(&xL, currentPositionX);
	_mm_storeh_pd(&xH, currentPositionX);
	_mm_storel_pd(&yL, currentPositionY);
	_mm_storeh_pd(&yH, currentPositionY);

	//quadrants I and II : quadrants III and IV
	const double thetaL = (yL > 0.0) ? atan2(yL, xL) : 2.0*M_PI - atan2(-yL, xL);
	const double thetaH = (yH > 0.0) ? atan2(yH, xH) : 2.0*M_PI - atan2(-yH, xH);

	const __m128d cosV = _mm_set_pd(cos(thetaH), cos(thetaL));
	const __m128d sinV = _mm_set_pd(sin(thetaH), sin(thetaL));

	_mm_store_pd(pFx, _mm_load_pd(pFx) - decayV*expR*cosV);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - decayV*expR*sinV);
}

void ConfinementForceVoid::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the confinement force void is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ConfinementForceVoid bit:
		forceFlags |= ConfinementForceVoidFlag;	// compound bitwise OR

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
                           6, const_cast<char *> ("[V/m^2] (ConfinementForceVoid)"), error);
	if (!*error)
		fits_write_key_dbl(file, const_cast<char *> ("decay"), decay, 
                           6, const_cast<char *> ("[1/m] (ConfinementForceVoid)"), error);
}

void ConfinementForceVoid::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConst"), &confine, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("decay"), &decay, NULL, error);
}
