/*===- ConfinementForceVoid.cpp - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ConfinementForceVoid.h"
	
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

	double rL, rH;
	_mm_storel_pd(&rL, r);
	_mm_storeh_pd(&rH, r);
	const __m128d expR = _mm_set_pd(exp(-decay*rH), exp(-decay*rL))/r;

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) - decayV*expR*currentPositionX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - decayV*expR*currentPositionY);
}

void ConfinementForceVoid::writeForce(fitsfile * const file, int * const error) const {
	ConfinementForce::writeForce(file, error);
	
	// add flag indicating that the confinement force void is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ConfinementForceVoid bit:
		forceFlags |= ConfinementForceFlag;
		forceFlags |= ConfinementForceVoidFlag;

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConst"), confine, 
                           6, const_cast<char *> ("[V/m^2] (ConfinementForceVoid)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("decay"), decay, 
                           6, const_cast<char *> ("[m^-1] (ConfinementForceVoid)"), error);
	}
}

void ConfinementForceVoid::readForce(fitsfile * const file, int * const error) {
	ConfinementForce::readForce(file, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConst"), &confine, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("decay"), &decay, NULL, error);
}
