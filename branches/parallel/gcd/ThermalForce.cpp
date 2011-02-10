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
#include <ctime>
#include "VectorCompatibility.h"

ThermalForce::ThermalForce(Cloud * const myCloud, const double redFactor) : Force(myCloud), mt(time(NULL)), 
evenRandCache(new RandCache[myCloud->n/2]), oddRandCache(new RandCache[myCloud->n/2]), 
evenRandGroup(dispatch_group_create()), oddRandGroup(dispatch_group_create()),
randQueue(dispatch_queue_create("com.DEMON.ThermalForce", NULL)), 
heatVal(redFactor) 
{
	dispatch_group_async(oddRandGroup, randQueue, ^{
		for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
			oddRandCache[i] = RandCache(_mm_set_pd(mt(), mt()), mt(), mt());
	});
}

ThermalForce::~ThermalForce()
{
	delete[] evenRandCache;
	delete[] oddRandCache;
	dispatch_release(evenRandGroup);
	dispatch_release(oddRandGroup);
	dispatch_release(randQueue);
}

void ThermalForce::force1(const double currentTime)
{
	dispatch_group_async(evenRandGroup, randQueue, ^{
		for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
			evenRandCache[i] = RandCache(_mm_set_pd(mt(), mt()), mt(), mt());
	});
	
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
	dispatch_apply(cloud->n/2, queue, ^(cloud_index currentParticle) {
		force(currentParticle*2, oddRandCache[currentParticle]);
	});
}

void ThermalForce::force2(const double currentTime)
{
	dispatch_group_async(oddRandGroup, randQueue, ^{
		for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
			oddRandCache[i] = RandCache(_mm_set_pd(mt(), mt()), mt(), mt());
	});
	
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
	dispatch_apply(cloud->n/2, queue, ^(cloud_index currentParticle) {
		force(currentParticle*2, evenRandCache[currentParticle]);
	});
}

void ThermalForce::force3(const double currentTime)
{
	dispatch_group_async(evenRandGroup, randQueue, ^{
		for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
			evenRandCache[i] = RandCache(_mm_set_pd(mt(), mt()), mt(), mt());
	});
	
	dispatch_group_wait(oddRandGroup, DISPATCH_TIME_FOREVER);
	dispatch_apply(cloud->n/2, queue, ^(cloud_index currentParticle) {
		force(currentParticle*2, oddRandCache[currentParticle]);
	});
}

void ThermalForce::force4(const double currentTime)
{
	dispatch_group_async(oddRandGroup, randQueue, ^{
		for (cloud_index i = 0, e = cloud->n/2; i < e; i++)
			oddRandCache[i] = RandCache(_mm_set_pd(mt(), mt()), mt(), mt());
	});
	
	dispatch_group_wait(evenRandGroup, DISPATCH_TIME_FOREVER);
	dispatch_apply(cloud->n/2, queue, ^(cloud_index currentParticle) {
		force(currentParticle*2, evenRandCache[currentParticle]);
	});
}

inline void ThermalForce::force(const cloud_index currentParticle, const RandCache &rc)
{	
	// MT random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*rc.r;
	const double thetaL = rc.l*2.0*M_PI;
	const double thetaH = rc.h*2.0*M_PI;
	
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	
	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(cos(thetaH), cos(thetaL)));
}

void ThermalForce::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the thermal force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ThermalForce bit:
		forceFlags |= ThermalForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
	{
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue"), heatVal, 6, const_cast<char *> ("[N] (ThermalForce)"), error);
	}
}

void ThermalForce::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
	{
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue"), &heatVal, NULL, error);
	}
}
