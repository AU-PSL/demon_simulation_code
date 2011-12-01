/*===- ThermalForce.cpp - libSimulation -=======================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "ThermalForce.h"
#include <cmath>

ThermalForce::ThermalForce(Cloud * const C, const double redFactor) 
: Force(C), evenRandCache(new RandCache[C->n/2]), oddRandCache(new RandCache[C->n/2]),
#ifdef DISPATCH_QUEUES
evenRandGroup(dispatch_group_create()), oddRandGroup(dispatch_group_create()),
randQueue(dispatch_queue_create("com.DEMON.ThermalForce", NULL)),
#endif
heatVal(redFactor) {
#ifdef DISPATCH_QUEUES
    dispatch_group_async(oddRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
        oddRandCache[i] = RandCache(_mm_set_pd(cloud->rands.uniformZeroToOne(), 
											   cloud->rands.uniformZeroToOne()), 
									cloud->rands.uniformZeroToTwoPi(), 
									cloud->rands.uniformZeroToTwoPi());
#ifdef DISPATCH_QUEUES
    });
#endif
}

ThermalForce::~ThermalForce() {
    delete[] evenRandCache;
    delete[] oddRandCache;

#ifdef DISPATCH_QUEUES
    dispatch_release(evenRandGroup);
	dispatch_release(oddRandGroup);
	dispatch_release(randQueue);
#endif
}

void ThermalForce::force1(const double currentTime) {
    (void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(evenRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
        evenRandCache[i] = RandCache(_mm_set_pd(cloud->rands.uniformZeroToOne(), 
												cloud->rands.uniformZeroToOne()), 
									 cloud->rands.uniformZeroToTwoPi(), 
									 cloud->rands.uniformZeroToTwoPi());
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, oddRandCache[currentParticle/2]);
    END_PARALLEL_FOR
}

void ThermalForce::force2(const double currentTime) {
    (void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(oddRandGroup, randQueue, ^{
#endif
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
        oddRandCache[i] = RandCache(_mm_set_pd(cloud->rands.uniformZeroToOne(), 
											   cloud->rands.uniformZeroToOne()), 
									cloud->rands.uniformZeroToTwoPi(), 
									cloud->rands.uniformZeroToTwoPi());
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
#endif

    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
        force(currentParticle, evenRandCache[currentParticle/2]);
    END_PARALLEL_FOR
}

void ThermalForce::force3(const double currentTime) {
    (void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(evenRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
        evenRandCache[i] = RandCache(_mm_set_pd(cloud->rands.uniformZeroToOne(), 
												cloud->rands.uniformZeroToOne()), 
									 cloud->rands.uniformZeroToTwoPi(), 
									 cloud->rands.uniformZeroToTwoPi());
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, oddRandCache[currentParticle/2]);
    END_PARALLEL_FOR
}

void ThermalForce::force4(const double currentTime) {
    (void)currentTime;
#ifdef DISPATCH_QUEUES
    dispatch_group_async(oddRandGroup, randQueue, ^{
#endif
    for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
        oddRandCache[i] = RandCache(_mm_set_pd(cloud->rands.uniformZeroToOne(), 
											   cloud->rands.uniformZeroToOne()), 
									cloud->rands.uniformZeroToTwoPi(), 
									cloud->rands.uniformZeroToTwoPi());
#ifdef DISPATCH_QUEUES
    });
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
#endif
    
    BEGIN_PARALLEL_FOR(currentParticle, numParticles, cloud->n, 2, static) 
		force(currentParticle, evenRandCache[currentParticle/2]);
    END_PARALLEL_FOR
}

// F = c*L : L is a uniformly distributed random number between 0 - 1 in a 
// random direction.
inline void ThermalForce::force(const cloud_index currentParticle, const RandCache &rc) {
	const doubleV thermV = _mm_set1_pd(heatVal)*rc.r;
	
	plusEqual_pd(cloud->forceX + currentParticle, thermV*_mm_set_pd(sin(rc.h), sin(rc.l))); // _mm_set_pd() is backwards
	plusEqual_pd(cloud->forceY + currentParticle, thermV*_mm_set_pd(sin(rc.h), sin(rc.l)));
}

void ThermalForce::writeForce(fitsfile * const file, int * const error) const {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the thermal force is used:
	if (!*error) {
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ThermalForce bit:
		forceFlags |= ThermalForceFlag;

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, 
                            const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue"), heatVal, 
                           6, const_cast<char *> ("[N] (ThermalForce)"), error);
}

void ThermalForce::readForce(fitsfile * const file, int * const error) {
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue"), &heatVal, NULL, error);
}
