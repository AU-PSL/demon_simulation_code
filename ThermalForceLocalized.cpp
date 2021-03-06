/**
* @file  ThermalForceLocalized.cpp
* @class ThermalForceLocalized ThermalForceLocalized.h
*
* @brief Computes a random force to model thermal effects
*		 where the force changes at a given radius
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "ThermalForceLocalized.h"
#include <cmath>

ThermalForceLocalized::ThermalForceLocalized(Cloud * const C, const double thermRed1, 
                                             const double thermRed2, const double specifiedRadius) 
: Force(C), heatingRadius(specifiedRadius), heatVal1(thermRed1), heatVal2(thermRed2), 
evenRandCache(new RandCache[C->n/DOUBLE_STRIDE]), oddRandCache(new RandCache[C->n/DOUBLE_STRIDE])
#ifdef DISPATCH_QUEUES
, evenRandGroup(dispatch_group_create()), oddRandGroup(dispatch_group_create()),
randQueue(dispatch_queue_create("com.DEMON.ThermalForceLocalized", NULL))
#endif
{
    for (cloud_index i = 0, e = cloud->n/DOUBLE_STRIDE; i < e; i++)
        oddRandCache[i] = RandCache(cloud->rands);
}

ThermalForceLocalized::~ThermalForceLocalized() {
    delete[] evenRandCache;
	delete[] oddRandCache;
    
#ifdef DISPATCH_QUEUES
	dispatch_release(evenRandGroup);
	dispatch_release(oddRandGroup);
	dispatch_release(randQueue);
#endif
}

void ThermalForceLocalized::force1(const double currentTime) {
    (void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(evenRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/DOUBLE_STRIDE; i < e; i++)
        evenRandCache[i] = RandCache(cloud->rands);
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
#endif   

	BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), 
              oddRandCache[currentParticle/DOUBLE_STRIDE]);
    END_PARALLEL_FOR
}

void ThermalForceLocalized::force2(const double currentTime) {
	(void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(oddRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/DOUBLE_STRIDE; i < e; i++)
        oddRandCache[i] = RandCache(cloud->rands);
#ifdef DISPATCH_QUEUES
	});
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), 
              evenRandCache[currentParticle/DOUBLE_STRIDE]);
    END_PARALLEL_FOR
}

void ThermalForceLocalized::force3(const double currentTime) {
	(void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(evenRandGroup, randQueue, ^{
#endif
        for (cloud_index i = 0, e = cloud->n/DOUBLE_STRIDE; i < e; i++)
            evenRandCache[i] = RandCache(cloud->rands);
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
    force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), 
          oddRandCache[currentParticle/DOUBLE_STRIDE]);
    END_PARALLEL_FOR
}

void ThermalForceLocalized::force4(const double currentTime) {
	(void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(oddRandGroup, randQueue, ^{
#endif
        for (cloud_index i = 0, e = cloud->n/DOUBLE_STRIDE; i < e; i++)
            oddRandCache[i] = RandCache(cloud->rands);
#ifdef DISPATCH_QUEUES
	});
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, DOUBLE_STRIDE, static) 
    force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), 
          evenRandCache[currentParticle/DOUBLE_STRIDE]);
    END_PARALLEL_FOR
}

// F = c1*L : if r > h_r 
// F = c2*L : if r < h_r 
// L is a uniformly distributed random number between 0 - 1 in a random 
// direction.
/**
* @brief Computes a thermal force with form F = c1*L : if r > h_r 
*										and F = c2*L : if r < h_r
*		 where L is a uniformly distributed random number between 0 - 1 in a 
*        random direction.
*
* @param[in] currentParticle The particle whose force is being computed
* @param[in] RC              RandCache struct from RandomNumbers.h
**/
inline void ThermalForceLocalized::force(const cloud_index currentParticle, const doubleV displacementX, 
                                         const doubleV displacementY, const RandCache &RC) {
	const doubleV radiusV = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	
	const int mask = movemask_pd(cmplt_pd(radiusV, heatingRadius));
    const doubleV thermV = mul_pd(select_pd(mask, heatVal1, heatVal2), RC.r);
	
	plusEqual_pd(cloud->forceX + currentParticle, thermV*randomCos(RC)); // _mm_set_pd() is backwards
	plusEqual_pd(cloud->forceY + currentParticle, thermV*randomSin(RC));
}

inline const doubleV ThermalForceLocalized::randomCos(const RandCache &RC) {
#ifdef __AVX__
    return _mm256_set_pd(cos(RC.r4), cos(RC.r3), cos(RC.r2), cos(RC.r1));
#else
    return _mm_set_pd(cos(RC.r2), cos(RC.r1));
#endif
}

inline const doubleV ThermalForceLocalized::randomSin(const RandCache &RC) {
#ifdef __AVX__
    return _mm256_set_pd(sin(RC.r4), sin(RC.r3), sin(RC.r2), sin(RC.r1));
#else
    return _mm_set_pd(sin(RC.r2), sin(RC.r1));
#endif 
}

void ThermalForceLocalized::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the localized thermal force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ThermalForce bit:
		forceFlags |= ThermalForceLocalizedFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error) {
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue1"), heatVal1, 
                           6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue2"), heatVal2, 
                           6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingRadius"), heatingRadius, 
                           6, const_cast<char *> ("[m] (ThermalForceLocalized)"), error);
	}
}

void ThermalForceLocalized::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error) {
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue1"), &heatVal1, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue2"), &heatVal2, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingRadius"), &heatingRadius, NULL, error);
	}
}
