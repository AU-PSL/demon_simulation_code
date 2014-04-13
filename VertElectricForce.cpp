/*===- VertElectricForce.cpp - libSimulation -=====================================
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

#include "VertElectricForce.h"

void VertElectricForce::force1(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->gety1_pd(currentParticle));
    END_PARALLEL_FOR
}

void VertElectricForce::force2(const double currentTime) {
    (void)currentTime;
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
        force(currentParticle, cloud->gety2_pd(currentParticle));
    END_PARALLEL_FOR
}

void VertElectricForce::force3(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->gety3_pd(currentParticle));
    END_PARALLEL_FOR
}

void VertElectricForce::force4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->gety4_pd(currentParticle));
    END_PARALLEL_FOR
}

// F = q*E*(y-top)/(bottom-top)
inline void VertElectricForce::force(const cloud_index currentParticle, const doubleV currentPositionY) {

	const doubleV cV = mul_pd(load_pd(cloud->charge + currentParticle), vertElectric);

        //const doubleV num = sub_pd(currentPositionY, top);
        //const double denom = top-bottom;
        //const doubleV frac = div_pd(num,denom);
	//const doubleV eV = mul_pd(cV,frac);

	const doubleV eV = mul_pd(cV, exp_pd(div_pd(currentPositionY,vertDec)));
	
	plusEqual_pd(cloud->forceY + currentParticle, eV);
}

void VertElectricForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the vertElectric force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add VertElectricForce bit:
		forceFlags |= VertElectricForceFlag;
		
		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("VertElectricFieldStrength"), vertElectric, 
                           6, const_cast<char *> ("[V/m^2] (VertElectricForce)"), error);
        	fits_write_key_dbl(file, const_cast<char *> ("verticalDecay"), vertDec, 
                           6, const_cast<char *> ("[m] (VertElectricForce)"), error);

}

void VertElectricForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("vertElectricFieldStrength"), &vertElectric, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("verticalDecay"), &vertDec, NULL, error);

}
