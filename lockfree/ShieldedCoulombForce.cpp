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
: Force(C), shielding(shieldingConstant) {}
ShieldedCoulombForce::~ShieldedCoulombForce() {}

void ShieldedCoulombForce::force1(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force1 does not fully support AVX."
#endif
    (void)currentTime;
    coulombForce1(0, cloud->n);
}

void ShieldedCoulombForce::coulombForce1(const cloud_index i, const cloud_index j) {
    const cloud_index n = j - i;
    if (n > cloud->n/16) {
        const cloud_index half = i + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce1(i, half);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce1(half, j);});
        block1(i, half, half, j);
    } else {
        for (cloud_index ci = i, e = j - 1; ci < e; ci += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx1_pd(i);
            const doubleV vy1 = cloud->gety1_pd(i);
            const doubleV vq1 = _mm_load_pd(cloud->charge + i);
            double x1, x2, y1, y2, q1, q2;
            _mm_storel_pd(&x1, vx1);
            _mm_storeh_pd(&x2, vx1);
            _mm_storel_pd(&y1, vy1);
            _mm_storeh_pd(&y2, vy1);
            _mm_storel_pd(&q1, vq1);
            _mm_storeh_pd(&q2, vq1);
        
            force(i, q1, q2, x1 - x2, y1 - y2);
            
            for (cloud_index cj = i + DOUBLE_STRIDE; cj < j; cj += DOUBLE_STRIDE) {
                double * const c = cloud->charge + cj;
                force(i, i, vq1, _mm_load_pd(c), vx1 - cloud->getx1_pd(cj), vy1 - cloud->gety1_pd(cj));
                forcer(i, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx1r_pd(cj), vy1 - cloud->gety1r_pd(cj));
            }
        }
    }
}

void ShieldedCoulombForce::block1(const cloud_index ilow, const cloud_index ihigh, 
                                  const cloud_index jlow, const cloud_index jhigh) {
    const cloud_index m = ihigh - ilow;
    const cloud_index n = jhigh - jlow;
    
    if (m*n > cloud->n*cloud->n/256) {
        const cloud_index mhalf = ilow + m/2;
        const cloud_index nhalf = jlow + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block1(mhalf, ihigh, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block1(ilow, mhalf, nhalf, jhigh);});
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block1(ilow, mhalf, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block1(mhalf, ihigh, nhalf, jhigh);});
    } else
        for (cloud_index i = ilow, ie = ihigh - 1; i < ie; i += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx1_pd(i);
            const doubleV vy1 = cloud->gety1_pd(i);
            
            const doubleV vq1 = _mm_load_pd(cloud->charge + ilow);
            
            for (cloud_index j = jlow, je = jhigh - 1; j < je; j += DOUBLE_STRIDE) {
                double * const c = cloud->charge + jlow;
                force(i, j, vq1, load_pd(c), vx1 - cloud->getx1_pd(j), vy1 - cloud->gety1_pd(j));
                forcer(i, j, vq1, _mm_loadr_pd(c), vx1 - cloud->getx1r_pd(j), vy1 - cloud->gety1r_pd(j));
            }
        }
}

void ShieldedCoulombForce::force2(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force2 does not fully support AVX."
#endif
    (void)currentTime;
    coulombForce2(0, cloud->n);
}

void ShieldedCoulombForce::coulombForce2(const cloud_index i, const cloud_index j) {
    const cloud_index n = j - i;
    if (n > cloud->n/16) {
        const cloud_index half = i + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce2(i, half);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce2(half, j);});
        block2(i, half, half, j);
    } else
        for (cloud_index ci = i, e = j - 1; ci < e; ci += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx2_pd(i);
            const doubleV vy1 = cloud->gety2_pd(i);
            const doubleV vq1 = _mm_load_pd(cloud->charge + i);
            double x1, x2, y1, y2, q1, q2;
            _mm_storel_pd(&x1, vx1);
            _mm_storeh_pd(&x2, vx1);
            _mm_storel_pd(&y1, vy1);
            _mm_storeh_pd(&y2, vy1);
            _mm_storel_pd(&q1, vq1);
            _mm_storeh_pd(&q2, vq1);
            
            force(i, q1, q2, x1 - x2, y1 - y2);
            
            for (cloud_index cj = i + DOUBLE_STRIDE; cj < j; cj += DOUBLE_STRIDE) {
                double * const c = cloud->charge + cj;
                force(i, i, vq1, _mm_load_pd(c), vx1 - cloud->getx2_pd(cj), vy1 - cloud->gety2_pd(cj));
                forcer(i, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx2r_pd(cj), vy1 - cloud->gety2r_pd(cj));
            }
        }
}

void ShieldedCoulombForce::block2(const cloud_index ilow, const cloud_index ihigh, 
                                  const cloud_index jlow, const cloud_index jhigh) {
    const cloud_index m = ihigh - ilow;
    const cloud_index n = jhigh - jlow;
    
    if (m*n > cloud->n*cloud->n/256) {
        const cloud_index mhalf = ilow + m/2;
        const cloud_index nhalf = jlow + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block2(mhalf, ihigh, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block2(ilow, mhalf, nhalf, jhigh);});
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block2(ilow, mhalf, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block2(mhalf, ihigh, nhalf, jhigh);});
    } else
       for (cloud_index i = ilow, ie = ihigh - 1; i < ie; i += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx2_pd(i);
            const doubleV vy1 = cloud->gety2_pd(i);
        
            const doubleV vq1 = _mm_load_pd(cloud->charge + ilow);
            
            for (cloud_index j = jlow, je = jhigh - 1; j < je; j += DOUBLE_STRIDE) {
                double * const c = cloud->charge + jlow;
                force(i, j, vq1, load_pd(c), vx1 - cloud->getx2_pd(j), vy1 - cloud->gety2_pd(j));
                forcer(i, j, vq1, _mm_loadr_pd(c), vx1 - cloud->getx2r_pd(j), vy1 - cloud->gety2r_pd(j));
            }
        }
}

void ShieldedCoulombForce::force3(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force3 does not fully support AVX."
#endif
    (void)currentTime;
    coulombForce3(0, cloud->n);
}

void ShieldedCoulombForce::coulombForce3(const cloud_index i, const cloud_index j) {
    const cloud_index n = j - i;
    if (n > cloud->n/16) {
        const cloud_index half = i + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce3(i, half);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce3(half, j);});
        block3(i, half, half, j);
    } else {
        for (cloud_index ci = i, e = j - 1; ci < e; ci += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx3_pd(i);
            const doubleV vy1 = cloud->gety3_pd(i);
            const doubleV vq1 = _mm_load_pd(cloud->charge + i);
            double x1, x2, y1, y2, q1, q2;
            _mm_storel_pd(&x1, vx1);
            _mm_storeh_pd(&x2, vx1);
            _mm_storel_pd(&y1, vy1);
            _mm_storeh_pd(&y2, vy1);
            _mm_storel_pd(&q1, vq1);
            _mm_storeh_pd(&q2, vq1);
            
            force(i, q1, q2, x1 - x2, y1 - y2);
            
            for (cloud_index cj = i + DOUBLE_STRIDE; cj < j; cj += DOUBLE_STRIDE) {
                double * const c = cloud->charge + cj;
                force(i, i, vq1, _mm_load_pd(c), vx1 - cloud->getx3_pd(cj), vy1 - cloud->gety3_pd(cj));
                forcer(i, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx3r_pd(cj), vy1 - cloud->gety3r_pd(cj));
            }
        }
    }
}

void ShieldedCoulombForce::block3(const cloud_index ilow, const cloud_index ihigh, 
                                  const cloud_index jlow, const cloud_index jhigh) {
    const cloud_index m = ihigh - ilow;
    const cloud_index n = jhigh - jlow;
    
    if (m*n > cloud->n*cloud->n/256) {
        const cloud_index mhalf = ilow + m/2;
        const cloud_index nhalf = jlow + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block3(mhalf, ihigh, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block3(ilow, mhalf, nhalf, jhigh);});
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block3(ilow, mhalf, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block3(mhalf, ihigh, nhalf, jhigh);});
    } else
        for (cloud_index i = ilow, ie = ihigh - 1; i < ie; i += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx3_pd(i);
            const doubleV vy1 = cloud->gety3_pd(i);
            
            const doubleV vq1 = _mm_load_pd(cloud->charge + ilow);
            
            for (cloud_index j = jlow, je = jhigh - 1; j < je; j += DOUBLE_STRIDE) {
                double * const c = cloud->charge + jlow;
                force(i, j, vq1, load_pd(c), vx1 - cloud->getx3_pd(j), vy1 - cloud->gety3_pd(j));
                forcer(i, j, vq1, _mm_loadr_pd(c), vx1 - cloud->getx3r_pd(j), vy1 - cloud->gety3r_pd(j));
            }
        }
}

void ShieldedCoulombForce::force4(const double currentTime) {
#ifdef __AVX__
#error "ShieldedCoulombForce::force4 does not fully support AVX."
#endif
    (void)currentTime;
    coulombForce4(0, cloud->n);
}

void ShieldedCoulombForce::coulombForce4(const cloud_index i, const cloud_index j) {
    const cloud_index n = j - i;
    if (n > cloud->n/16) {
        const cloud_index half = i + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce4(i, half);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{coulombForce4(half, j);});
        block4(i, half, half, j);
    } else {
        for (cloud_index ci = i, e = j - 1; ci < e; ci += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx4_pd(i);
            const doubleV vy1 = cloud->gety4_pd(i);
            const doubleV vq1 = _mm_load_pd(cloud->charge + i);
            double x1, x2, y1, y2, q1, q2;
            _mm_storel_pd(&x1, vx1);
            _mm_storeh_pd(&x2, vx1);
            _mm_storel_pd(&y1, vy1);
            _mm_storeh_pd(&y2, vy1);
            _mm_storel_pd(&q1, vq1);
            _mm_storeh_pd(&q2, vq1);
            
            force(i, q1, q2, x1 - x2, y1 - y2);
            
            for (cloud_index cj = i + DOUBLE_STRIDE; cj < j; cj += DOUBLE_STRIDE) {
                double * const c = cloud->charge + cj;
                force(i, i, vq1, _mm_load_pd(c), vx1 - cloud->getx4_pd(cj), vy1 - cloud->gety4_pd(cj));
                forcer(i, i, vq1, _mm_loadr_pd(c), vx1 - cloud->getx4r_pd(cj), vy1 - cloud->gety4r_pd(cj));
            }
        }
    }
}

void ShieldedCoulombForce::block4(const cloud_index ilow, const cloud_index ihigh, 
                                  const cloud_index jlow, const cloud_index jhigh) {
    const cloud_index m = ihigh - ilow;
    const cloud_index n = jhigh - jlow;
    
    if (m*n > cloud->n*cloud->n/256) {
        const cloud_index mhalf = ilow + m/2;
        const cloud_index nhalf = jlow + n/2;
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block4(mhalf, ihigh, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block4(ilow, mhalf, nhalf, jhigh);});
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block4(ilow, mhalf, jlow, nhalf);});
        dispatch_sync(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{block4(mhalf, ihigh, nhalf, jhigh);});
    } else
        for (cloud_index i = ilow, ie = ihigh - 1; i < ie; i += DOUBLE_STRIDE) {
            const doubleV vx1 = cloud->getx4_pd(i);
            const doubleV vy1 = cloud->gety4_pd(i);
            
            const doubleV vq1 = _mm_load_pd(cloud->charge + ilow);
            
            for (cloud_index j = jlow, je = jhigh - 1; j < je; j += DOUBLE_STRIDE) {
                double * const c = cloud->charge + jlow;
                force(i, j, vq1, load_pd(c), vx1 - cloud->getx4_pd(j), vy1 - cloud->gety4_pd(j));
                forcer(i, j, vq1, _mm_loadr_pd(c), vx1 - cloud->getx4r_pd(j), vy1 - cloud->gety4r_pd(j));
            }
        }
}

// F_i,i+1 = e0*q_i*q_i+1/(|r_i - r_i+1|)^2*Exp(-s*|r_i - r_i+1|)*(1 + c*|r_i - r_i+1|)
// Calculates inteaction between i and i+1 particles.
inline void ShieldedCoulombForce::force(const cloud_index currentParticle,
                                        const double iCharge, const double jCharge,
                                        const double displacementX, const double displacementY) {
	// Calculate displacement between particles.
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement*shielding;
    
	if (valExp < 10.0) {// restrict to 10*(ion debye length)
		// calculate force
		const double forceC = coulomb*iCharge*jCharge*(1.0 + valExp)
            /(displacement*displacement*displacement*exp(valExp));
        
        const double forceX = forceC*displacementX;
        const double forceY = forceC*displacementY;
		plusEqual_pd(cloud->forceX + currentParticle, _mm_set_pd(-forceX, forceX));
        plusEqual_pd(cloud->forceY + currentParticle, _mm_set_pd(-forceY, forceY));
	}
}

// F_i,j = e0*q_i*q_j/(|r_i - r_j|)^2*Exp(-s*|r_i - r_j|)*(1 + c*|r_i - r_j|)
// Calculates inteaction between (i,i+1) and (j,j+1) particles.
inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const doubleV currentCharge, const doubleV iCharge,
                                        const doubleV displacementX, const doubleV displacementY) {
	// Calculate displacement between particles.
	const doubleV displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const doubleV valExp = displacement*_mm_set1_pd(shielding);
	
	const int mask = _mm_movemask_pd(_mm_cmplt_pd(valExp, _mm_set1_pd(10.0)));
	if (!mask)
		return;
	
	double expl = 0.0, exph = 0.0;
    if (mask & 1) {
        double expVal;
        _mm_storel_pd(&expVal, valExp);
        expl = exp(-expVal);
    }
    if (mask & 2) {
        double expVal;
        _mm_storeh_pd(&expVal, valExp);
        exph = exp(-expVal);
    }
    
    // calculate force
	const doubleV forceC = _mm_set1_pd(coulomb)*currentCharge*iCharge*(_mm_set1_pd(1.0) + valExp)*_mm_set_pd(exph, expl)
    /(displacement*displacement*displacement);
    const doubleV forcevX = forceC*displacementX;
	const doubleV forcevY = forceC*displacementY;
    
	plusEqual_pd(cloud->forceX + currentParticle, forcevX);
	plusEqual_pd(cloud->forceY + currentParticle, forcevY);

	// equal and opposite force:
    minusEqual_pd(cloud->forceX + iParticle, forcevX);
	minusEqual_pd(cloud->forceY + iParticle, forcevY);
}

// F_i,j = e0*q_i*q_j/(|r_i - r_j|)^2*Exp(-s*|r_i - r_j|)*(1 + c*|r_i - r_j|)
// Calculates inteaction between (i,i+1) and (j+1,j) particles.
inline void ShieldedCoulombForce::forcer(const cloud_index currentParticle, const cloud_index iParticle, 
                                         const doubleV currentCharge, const doubleV iCharge,
                                         const doubleV displacementX, const doubleV displacementY) {
	// Calculate displacement between particles.
	const doubleV displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const doubleV valExp = displacement*_mm_set1_pd(shielding);
	
	const int mask = _mm_movemask_pd(_mm_cmplt_pd(valExp, _mm_set1_pd(10.0)));
	if (!mask)
		return;
    
    double expl = 0.0, exph = 0.0;
    if (mask & 1) {
        double expVal;
        _mm_storel_pd(&expVal, valExp);
        expl = exp(-expVal);
    }
    if (mask & 2) {
        double expVal;
        _mm_storeh_pd(&expVal, valExp);
        exph = exp(-expVal);
    }
    
    // calculate force
	const doubleV forceC = _mm_set1_pd(coulomb)*currentCharge*iCharge*(_mm_set1_pd(1.0) + valExp)*_mm_set_pd(exph, expl)
    /(displacement*displacement*displacement);
	const doubleV forcevX = forceC*displacementX;
	const doubleV forcevY = forceC*displacementY;
    
    plusEqual_pd(cloud->forceX + currentParticle, forcevX);
	plusEqual_pd(cloud->forceY + currentParticle, forcevY);

	// equal and opposite force:
    minusEqualr_pd(cloud->forceX + iParticle, forcevX);
	minusEqualr_pd(cloud->forceY + iParticle, forcevY);
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

inline void ShieldedCoulombForce::plusEqualr_pd(double * const a, const doubleV b) {
	_mm_storer_pd(a, _mm_add_pd(_mm_loadr_pd(a), b));
}

inline void ShieldedCoulombForce::minusEqualr_pd(double * const a, const doubleV b) {
	_mm_storer_pd(a, _mm_sub_pd(_mm_loadr_pd(a), b));
}
