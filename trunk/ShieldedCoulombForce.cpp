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

#ifdef DISPATCH_QUEUES
#define LOOP_END(num) num
#else
#define LOOP_END(num) num - 1
#endif

const double ShieldedCoulombForce::coulomb = 1.0/(4.0*M_PI*Cloud::epsilon0);

ShieldedCoulombForce::ShieldedCoulombForce(Cloud * const C, const double shieldingConstant)
: Force(C), shielding(shieldingConstant) SEMAPHORES_MALLOC(C->n/DOUBLE_STRIDE) {
    SEMAPHORES_INIT(cloud->n/DOUBLE_STRIDE)
}
ShieldedCoulombForce::~ShieldedCoulombForce() {
    SEMAPHORES_FREE(cloud->n/DOUBLE_STRIDE)
}

// FIXME: When changing this over to AVX to simplify, change force methods to
// triangle and block forces. Triangle forces cover the outter loop force. Block
// covers the inner loop forces. AVX specific differences would go into that 
// section.
void ShieldedCoulombForce::force1(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force1 does not fully support AVX."
#endif
    (void)currentTime;
    const cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, LOOP_END(numParticles), DOUBLE_STRIDE, dynamic)
        const doubleV vx1 = cloud->getx1_pd(currentParticle);
        const doubleV vy1 = cloud->gety1_pd(currentParticle);
        const doubleV vq1 = load_pd(cloud->charge + currentParticle);
        double x1, x2, y1, y2, q1, q2;
        _mm_storel_pd(&x1, vx1);
        _mm_storeh_pd(&x2, vx1);
        _mm_storel_pd(&y1, vy1);
        _mm_storeh_pd(&y2, vy1);
        _mm_storel_pd(&q1, vq1);
        _mm_storeh_pd(&q2, vq1);
                 
        force(currentParticle, q1, q2, x1 - x2, y1 - y2);
        for (cloud_index i = currentParticle + DOUBLE_STRIDE; i < numParticles; i += DOUBLE_STRIDE) {
			double * const c = cloud->charge + i;
            force(currentParticle, i, vq1, load_pd(c), sub_pd(vx1, cloud->getx1_pd(i)), sub_pd(vy1, cloud->gety1_pd(i)));
            forcer(currentParticle, i, vq1, _mm_loadr_pd(c), sub_pd(vx1, cloud->getx1r_pd(i)), sub_pd(vy1, cloud->gety1r_pd(i)));
        }
    END_PARALLEL_FOR
}

void ShieldedCoulombForce::force2(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force2 does not fully support AVX."
#endif
    (void)currentTime;
	const cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, LOOP_END(numParticles), DOUBLE_STRIDE, dynamic)
		const doubleV vx1 = cloud->getx2_pd(currentParticle);
		const doubleV vy1 = cloud->gety2_pd(currentParticle);
		const doubleV vq1 = load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + DOUBLE_STRIDE; i < numParticles; i += DOUBLE_STRIDE) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, load_pd(c), sub_pd(vx1, cloud->getx2_pd(i)), sub_pd(vy1, cloud->gety2_pd(i)));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), sub_pd(vx1, cloud->getx2r_pd(i)), sub_pd(vy1, cloud->gety2r_pd(i)));
		}
	END_PARALLEL_FOR
}

void ShieldedCoulombForce::force3(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force3 does not fully support AVX."
#endif
    (void)currentTime;
    const cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, LOOP_END(numParticles), DOUBLE_STRIDE, dynamic)
		const doubleV vx1 = cloud->getx3_pd(currentParticle);
		const doubleV vy1 = cloud->gety3_pd(currentParticle);
		const doubleV vq1 = load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + DOUBLE_STRIDE; i < numParticles; i += DOUBLE_STRIDE) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, load_pd(c), sub_pd(vx1, cloud->getx3_pd(i)), sub_pd(vy1, cloud->gety3_pd(i)));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), sub_pd(vx1, cloud->getx3r_pd(i)), sub_pd(vy1, cloud->gety3r_pd(i)));
		}
	END_PARALLEL_FOR
}

void ShieldedCoulombForce::force4(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force4 does not fully support AVX."
#endif
    (void)currentTime;
	const cloud_index numParticles = cloud->n;
    BEGIN_PARALLEL_FOR(currentParticle, e, LOOP_END(numParticles), DOUBLE_STRIDE, dynamic)
		const doubleV vx1 = cloud->getx4_pd(currentParticle);
		const doubleV vy1 = cloud->gety4_pd(currentParticle);
		const doubleV vq1 = load_pd(cloud->charge + currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + DOUBLE_STRIDE; i < numParticles; i += DOUBLE_STRIDE) {
			double * const c = cloud->charge + i;
			force(currentParticle, i, vq1, load_pd(c), sub_pd(vx1, cloud->getx4_pd(i)), sub_pd(vy1, cloud->gety4_pd(i)));
			forcer(currentParticle, i, vq1, _mm_loadr_pd(c), sub_pd(vx1, cloud->getx4r_pd(i)), sub_pd(vy1, cloud->gety4r_pd(i)));
		}
	END_PARALLEL_FOR
}

// F_i,i+1 = e0*q_i*q_i+1/(|r_i - r_i+1|)^2*Exp(-s*|r_i - r_i+1|)*(1 + c*|r_i - r_i+1|)
// Calculates inteaction between i and i+1 particles.
inline void ShieldedCoulombForce::force(const cloud_index currentParticle,
                                        const double currentCharge, const double iCharge,
                                        const double displacementX, const double displacementY) {
	// Calculate displacement between particles.
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement*shielding;

	if (valExp < 10.0) {// restrict to 10*(ion debye length)
		// calculate force
		const double forceC = coulomb*currentCharge*iCharge*(1.0 + valExp)
            /(displacement*displacement*displacement*exp(valExp));

        const double forceX = forceC*displacementX;
        const double forceY = forceC*displacementY;
        SEMAPHORE_WAIT(currentParticle/DOUBLE_STRIDE)
		plusEqual_pd(cloud->forceX + currentParticle, _mm_set_pd(-forceX, forceX));
        plusEqual_pd(cloud->forceY + currentParticle, _mm_set_pd(-forceY, forceY));
        SEMAPHORE_SIGNAL(currentParticle/DOUBLE_STRIDE)
	}
}

// F_i,j = e0*q_i*q_j/(|r_i - r_j|)^2*Exp(-s*|r_i - r_j|)*(1 + c*|r_i - r_j|)
// Calculates inteaction between (i,i+1) and (j,j+1) particles.
inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const doubleV currentCharge, const doubleV iCharge,
                                        const doubleV displacementX, const doubleV displacementY) {
	// Calculate displacement between particles.
	const doubleV displacement = sqrt_pd(add_pd(displacementX*displacementX, displacementY*displacementY));
	const doubleV valExp = displacement*set1_pd(shielding);
	
	const int mask = movemask_pd(_mm_cmplt_pd(valExp, set1_pd(10.0)));
	if (!mask)
		return;
    
    // calculate force
	const doubleV forceC = set1_pd(coulomb)*currentCharge*iCharge*(set1_pd(1.0) + valExp)*exp_pd(mask, valExp)
        /(displacement*displacement*displacement);
    const doubleV forcevX = forceC*displacementX;
	const doubleV forcevY = forceC*displacementY;

	SEMAPHORE_WAIT(currentParticle/DOUBLE_STRIDE)
	plusEqual_pd(cloud->forceX + currentParticle, forcevX);
	plusEqual_pd(cloud->forceY + currentParticle, forcevY);
    SEMAPHORE_SIGNAL(currentParticle/DOUBLE_STRIDE)

    SEMAPHORE_WAIT(iParticle/DOUBLE_STRIDE)
	// equal and opposite force:
    minusEqual_pd(cloud->forceX + iParticle, forcevX);
	minusEqual_pd(cloud->forceY + iParticle, forcevY);
    SEMAPHORE_SIGNAL(iParticle/DOUBLE_STRIDE)
}

// F_i,j = e0*q_i*q_j/(|r_i - r_j|)^2*Exp(-s*|r_i - r_j|)*(1 + c*|r_i - r_j|)
// Calculates inteaction between (i,i+1) and (j+1,j) particles.
inline void ShieldedCoulombForce::forcer(const cloud_index currentParticle, const cloud_index iParticle, 
                                         const doubleV currentCharge, const doubleV iCharge,
                                         const doubleV displacementX, const doubleV displacementY) {
	// Calculate displacement between particles.
	const doubleV displacement = sqrt_pd(add_pd(displacementX*displacementX, displacementY*displacementY));
	const doubleV valExp = displacement*set1_pd(shielding);
	
	const int mask = movemask_pd(_mm_cmplt_pd(valExp, set1_pd(10.0)));
	if (!mask)
		return;
    
    // calculate force
	const doubleV forceC = set1_pd(coulomb)*currentCharge*iCharge*(set1_pd(1.0) + valExp)*exp_pd(mask, valExp)
        /(displacement*displacement*displacement);
	const doubleV forcevX = forceC*displacementX;
	const doubleV forcevY = forceC*displacementY;

    SEMAPHORE_WAIT(currentParticle/DOUBLE_STRIDE)
    plusEqual_pd(cloud->forceX + currentParticle, forcevX);
	plusEqual_pd(cloud->forceY + currentParticle, forcevY);
    SEMAPHORE_SIGNAL(currentParticle/DOUBLE_STRIDE)
    
    SEMAPHORE_WAIT(iParticle/DOUBLE_STRIDE)
	// equal and opposite force:
    minusEqualr_pd(cloud->forceX + iParticle, forcevX);
	minusEqualr_pd(cloud->forceY + iParticle, forcevY);
    SEMAPHORE_SIGNAL(iParticle/DOUBLE_STRIDE)
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
		forceFlags |= ShieldedCoulombForceFlag;

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

inline doubleV exp_pd(const int mask, const doubleV a) {
	double expl = 0.0, exph = 0.0;
    if (mask & 1) {
        double expVal;
        _mm_storel_pd(&expVal, a);
        expl = exp(-expVal);
    }
    if (mask & 2) {
        double expVal;
        _mm_storeh_pd(&expVal, a);
        exph = exp(-expVal);
    }
	return _mm_set_pd(exph, expl);
}

inline void ShieldedCoulombForce::plusEqualr_pd(double * const a, const doubleV b) {
	_mm_storer_pd(a, _mm_add_pd(_mm_loadr_pd(a), b));
}

inline void ShieldedCoulombForce::minusEqualr_pd(double * const a, const doubleV b) {
	_mm_storer_pd(a, _mm_sub_pd(_mm_loadr_pd(a), b));
}
