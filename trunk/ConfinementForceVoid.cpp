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
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force2(const double currentTime) {
	ConfinementForce::force2(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force3(const double currentTime) {
	ConfinementForce::force3(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void ConfinementForceVoid::force4(const double currentTime) {
	ConfinementForce::force4(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

// F = F_c - d*Exp(-d*r)*r
inline void ConfinementForceVoid::force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY) {
	const doubleV decayV = mul_pd(load_pd(cloud->charge + currentParticle), -decay);
    const doubleV r = length_pd(currentPositionX, currentPositionY);
	const doubleV expR = div_pd(exp_pd(mul_pd(r, -decay)), r);
	
	plusEqual_pd(cloud->forceX + currentParticle, mul_pd(mul_pd(decayV, expR), currentPositionX));
	plusEqual_pd(cloud->forceY + currentParticle, mul_pd(mul_pd(decayV, expR), currentPositionY));
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

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("decay"), decay, 
                           6, const_cast<char *> ("[m^-1] (ConfinementForceVoid)"), error);
}

void ConfinementForceVoid::readForce(fitsfile * const file, int * const error) {
	ConfinementForce::readForce(file, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("decay"), &decay, NULL, error);
}
