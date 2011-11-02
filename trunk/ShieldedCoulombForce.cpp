/*===- ShieldedCoulombForce.cpp - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ShieldedCoulombForce.h"
#include <cmath>

const double ShieldedCoulombForce::coulomb = 1.0/(4.0*M_PI*Cloud::epsilon0);

ShieldedCoulombForce::ShieldedCoulombForce(Cloud * const C, const double shieldingConstant)
: Force(C), shielding(shieldingConstant) SEMAPHORES_MALLOC(C->n/2) {
    SEMAPHORES_INIT(cloud->n/2)
}
ShieldedCoulombForce::~ShieldedCoulombForce() {
    SEMAPHORES_FREE(cloud->n/2)
}

void ShieldedCoulombForce::force1(const double currentTime) {
    (void)currentTime;
    cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, numParticles - 1, 2, dynamic)
        const __m128d vx1 = cloud->getx1_pd(currentParticle);
        const __m128d vy1 = cloud->gety1_pd(currentParticle);
        const __m128d vq1 = _mm_load_pd(cloud->charge + currentParticle);
        double x1, x2, y1, y2, q1, q2;
        _mm_storel_pd(&x1, vx1);
        _mm_storeh_pd(&x2, vx1);
        _mm_storel_pd(&y1, vy1);
        _mm_storeh_pd(&y2, vy1);
        _mm_storel_pd(&q1, vq1);
        _mm_storeh_pd(&q2, vq1);
                 
        force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
                 
        for (cloud_index i = currentParticle + 2; i < numParticles; i += 2) {
			double * const c = cloud->charge + i;
            force(currentParticle, i, vq1, _mm_load_pd(c), vx1 - cloud->getx1_pd(i), vy1 - cloud->gety1_pd(i));
            forcer(currentParticle, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx1r_pd(i), vy1 - cloud->gety1r_pd(i));
        }
    END_PARALLEL_FOR
}

void ShieldedCoulombForce::force2(const double currentTime) {
    (void)currentTime;
	cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, numParticles - 1, 2, dynamic)
		const __m128d vx1 = cloud->getx2_pd(currentParticle);
		const __m128d vy1 = cloud->gety2_pd(currentParticle);
		const __m128d vq1 = _mm_load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, _mm_load_pd(c), vx1 - cloud->getx2_pd(i), vy1 - cloud->gety2_pd(i));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx2r_pd(i), vy1 - cloud->gety2r_pd(i));
		}
	END_PARALLEL_FOR
}

void ShieldedCoulombForce::force3(const double currentTime) {
    (void)currentTime;
    cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, numParticles - 1, 2, dynamic)
		const __m128d vx1 = cloud->getx3_pd(currentParticle);
		const __m128d vy1 = cloud->gety3_pd(currentParticle);
		const __m128d vq1 = _mm_load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, _mm_load_pd(c), vx1 - cloud->getx3_pd(i), vy1 - cloud->gety3_pd(i));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx3r_pd(i), vy1 - cloud->gety3r_pd(i));
		}
	END_PARALLEL_FOR
}

void ShieldedCoulombForce::force4(const double currentTime) {
    (void)currentTime;
	cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, numParticles - 1, 2, dynamic)
		const __m128d vx1 = cloud->getx4_pd(currentParticle);
		const __m128d vy1 = cloud->gety4_pd(currentParticle);
		const __m128d vq1 = _mm_load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, _mm_load_pd(c), vx1 - cloud->getx4_pd(i), vy1 - cloud->gety4_pd(i));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx4r_pd(i), vy1 - cloud->gety4r_pd(i));
		}
	END_PARALLEL_FOR
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const double currentCharge, const double iCharge,
                                        const double displacementX, const double displacementY) {
	// Calculate displacement between particles.
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement*shielding;

	if (valExp < 10.0) {// restrict to 10*(ion debye length)
		// calculate force
		const double coefficient = coulomb/(displacement*exp(valExp));
        const double forceC = currentCharge*iCharge*coefficient*(1.0 + valExp)/(displacement*displacement);
        
        SEMAPHORE_WAIT(currentParticle/2)
		cloud->forceX[currentParticle] += forceC*displacementX;
		cloud->forceY[currentParticle] += forceC*displacementY;

		// equal and opposite force:
		cloud->forceX[iParticle] -= forceC*displacementX;
		cloud->forceY[iParticle] -= forceC*displacementY;
        SEMAPHORE_SIGNAL(currentParticle/2)
	}
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const __m128d currentCharge, const __m128d iCharge,
                                        const __m128d displacementX, const __m128d displacementY) {
	// Calculate displacement between particles.
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement*_mm_set1_pd(shielding);
	
	const int mask = _mm_movemask_pd(_mm_cmplt_pd(valExp, _mm_set1_pd(10.0)));
	if (!mask)
		return;
	
	double valExpL = 0.0, valExpH = 0.0;
	if (mask & 1) _mm_storel_pd(&valExpL, valExp);
	if (mask & 2) _mm_storeh_pd(&valExpH, valExp);
	const __m128d expv = _mm_set_pd((mask & 2) ? exp(-valExpH) : 0.0, // _mm_set_pd is backwards
									(mask & 1) ? exp(-valExpL) : 0.0);

    // calculate force
	const __m128d coefficient = _mm_set1_pd(coulomb)/displacement*expv;
    const __m128d forceC = currentCharge*iCharge*coefficient*(_mm_set1_pd(1.0) + valExp)/(displacement*displacement);
    const __m128d forcevX = forceC*displacementX;
	const __m128d forcevY = forceC*displacementY;
    
    double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
    SEMAPHORE_WAIT(currentParticle/2)
    _mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);
    SEMAPHORE_SIGNAL(currentParticle/2)

    pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle;
    SEMAPHORE_WAIT(iParticle/2)
	// equal and opposite force:
    _mm_store_pd(pFx, _mm_load_pd(pFx) - forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - forcevY);
    SEMAPHORE_SIGNAL(iParticle/2)
}

inline void ShieldedCoulombForce::forcer(const cloud_index currentParticle, const cloud_index iParticle, 
                                         const __m128d currentCharge, const __m128d iCharge,
                                         const __m128d displacementX, const __m128d displacementY) {
	// Calculate displacement between particles.
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement*_mm_set1_pd(shielding);
	
	const int mask = _mm_movemask_pd(_mm_cmplt_pd(valExp, _mm_set1_pd(10.0)));
	if (!mask)
		return;
	
	double valExpL = 0.0, valExpH = 0.0;
	if (mask & 1) _mm_storel_pd(&valExpL, valExp);
	if (mask & 2) _mm_storeh_pd(&valExpH, valExp);
	const __m128d expv = _mm_set_pd((mask & 2) ? exp(-valExpH) : 0.0, // _mm_set_pd is backwards
									(mask & 1) ? exp(-valExpL) : 0.0);
	
    // calculate force
	const __m128d coefficient = _mm_set1_pd(coulomb)/displacement*expv;
	const __m128d forceC = currentCharge*iCharge*coefficient*(_mm_set1_pd(1.0) + valExp)/(displacement*displacement);
	const __m128d forcevX = forceC*displacementX;
	const __m128d forcevY = forceC*displacementY;
    
    double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
    SEMAPHORE_WAIT(currentParticle/2)
    _mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);
    SEMAPHORE_SIGNAL(currentParticle/2)
    
    pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle;
    SEMAPHORE_WAIT(iParticle/2)
	// equal and opposite force:
	_mm_storer_pd(pFx, _mm_loadr_pd(pFx) - forcevX);
	_mm_storer_pd(pFy, _mm_loadr_pd(pFy) - forcevY);
    SEMAPHORE_SIGNAL(iParticle/2)
}

void ShieldedCoulombForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the shielded Coulomb force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ShieldedCoulombForce bit:
		forceFlags |= ShieldedCoulombForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword.
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("shieldingConstant"), shielding, 
                           6, const_cast<char *> ("[m^-1] (ShieldedCoulombForce)"), error);
}

void ShieldedCoulombForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("shieldingConstant"), &shielding, NULL, error);
}
