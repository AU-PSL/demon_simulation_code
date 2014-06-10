/*===- ElectricForce.cpp - libSimulation -=====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ElectricForce.h"

void ElectricForce::force1(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void ElectricForce::force2(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void ElectricForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void ElectricForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

// F = q*E*(exp(-r/R))
inline void ElectricForce::force(const cloud_index currentParticle, const doubleV currentPositionX, 
                                    const doubleV currentPositionY) {
	const doubleV cV = mul_pd(load_pd(cloud->charge + currentParticle), electric);

        const doubleV PX = mul_pd(currentPositionX,currentPositionX);
        const doubleV PY = mul_pd(currentPositionY,currentPositionY);
        const doubleV Ra = add_pd(PX,PY);
        const doubleV R  = sqrt_pd(Ra);
        const doubleV rad= set1_pd(radius * (-1));
        const doubleV eV = mul_pd(cV, exp_pd(div_pd(R,rad)));
	
	plusEqual_pd(cloud->forceX + currentParticle, mul_pd(eV, currentPositionX));
	plusEqual_pd(cloud->forceY + currentParticle, mul_pd(eV, currentPositionY));

}

void ElectricForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the electric force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ElectricForce bit:
		forceFlags |= ElectricForceFlag;
		
		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("ElectricFieldStrength"), electric, 
                           6, const_cast<char *> ("[V/m^2] (ElectricForce)"), error);
        	fits_write_key_dbl(file, const_cast<char *> ("plasmaRadius"), radius, 
                           6, const_cast<char *> ("[m] (ElectricForce)"), error);

}

void ElectricForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("electricFieldStrength"), &electric, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("plasmaRadius"), &radius, NULL, error);

}
