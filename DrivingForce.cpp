/*===- DrivingForce.cpp - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "DrivingForce.h"
#include <cmath>

const double DrivingForce::waveNum = 2.0*M_PI/0.002; // wavelength = 2mm
const double DrivingForce::angFreq = 2.0*M_PI*10.0; // 10Hz

void DrivingForce::force1(const double currentTime) {
	const __m128d vtime = _mm_set1_pd(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, vtime, cloud->getx1_pd(currentParticle));
    END_PARALLEL_FOR
}

void DrivingForce::force2(const double currentTime) {
	const __m128d vtime = _mm_set1_pd(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, vtime, cloud->getx2_pd(currentParticle));
    END_PARALLEL_FOR
}

void DrivingForce::force3(const double currentTime) {
	const __m128d vtime = _mm_set1_pd(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, vtime, cloud->getx3_pd(currentParticle));
    END_PARALLEL_FOR
}

void DrivingForce::force4(const double currentTime) {
	const __m128d vtime = _mm_set1_pd(currentTime);
	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static)
		force(currentParticle, vtime, cloud->getx4_pd(currentParticle));
    END_PARALLEL_FOR
}

// F = a*sin(k*x - w*t)*Exp(-(x + x0)^2/b)
inline void DrivingForce::force(const cloud_index currentParticle, const __m128d currentTime, const __m128d currentPositionX) {
	// NOTE: This is different than the equation listed in the paper. The paper 
	// is incorrect.
	const __m128d distV = currentPositionX - _mm_set1_pd(shift);
	const __m128d sinArg = _mm_set1_pd(waveNum)*currentPositionX - _mm_set1_pd(angFreq)*currentTime;
	const __m128d expArg = -distV*distV/_mm_set1_pd(driveConst);
	
	// no SIMD trig instructions; break vectors and perform separately:
	double sinArgL, sinArgH, expArgL, expArgH;
	_mm_storel_pd(&sinArgL, sinArg);
	_mm_storeh_pd(&sinArgH, sinArg);
	_mm_storel_pd(&expArgL, expArg);
	_mm_storeh_pd(&expArgH, expArg);

	plusEqual_pd(cloud->forceX + currentParticle, 
				 _mm_set1_pd(amplitude)*_mm_set_pd(sin(sinArgH), sin(sinArgL))*_mm_set_pd(exp(expArgH), exp(expArgL))); // _mm_set_pd() is backwards
}

void DrivingForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the driving force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add DrivingForce bit:
		forceFlags |= DrivingForceFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("drivingAmplitude"), amplitude, 
                           6, const_cast<char *> ("[N] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingConst"), driveConst, 
                           6, const_cast<char *> ("[m^2] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingShift"), shift, 
                           6, const_cast<char *> ("[m] (DrivingForce)"), error);
	}
}

void DrivingForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error) {
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("drivingAmplitude"), &amplitude, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingConst"), &driveConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingShift"), &shift, NULL, error);
	}
}
