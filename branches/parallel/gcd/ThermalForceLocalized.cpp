/*===- ThermalForceLocalized.cpp - libSimulation -==============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "ThermalForceLocalized.h"
#include <ctime>
#include <cmath>

ThermalForceLocalized::ThermalForceLocalized(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius) 
: Force(myCloud), mt(time(NULL)), heatingRadius(specifiedRadius), heatVal1(thermRed1), heatVal2(thermRed2), semaphore(dispatch_semaphore_create(1)) {}

ThermalForceLocalized::~ThermalForceLocalized()
{
	dispatch_release(semaphore);
}

void ThermalForceLocalized::force1(const double currentTime)
{
	dispatch_apply(cloud->n/2, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t currentParticle) {
		currentParticle *= 2; 
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
	});
}

void ThermalForceLocalized::force2(const double currentTime)
{
	dispatch_apply(cloud->n/2, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t currentParticle) {
		currentParticle *= 2; 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
	});
}

void ThermalForceLocalized::force3(const double currentTime)
{
	dispatch_apply(cloud->n/2, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t currentParticle) {
		currentParticle *= 2; 
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
	});
}

void ThermalForceLocalized::force4(const double currentTime)
{
	dispatch_apply(cloud->n/2, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(size_t currentParticle) {
		currentParticle *= 2; 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
	});
}

inline void ThermalForceLocalized::force(const cloud_index currentParticle, const __m128d displacementX, const __m128d displacementY)
{
	const __m128d radiusV = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
	const double thetaL = mt()*2.0*M_PI;
	const double thetaH = mt()*2.0*M_PI;
	const __m128d randV = _mm_set_pd(mt(), mt());
	dispatch_semaphore_signal(semaphore);
	
	double rL, rH;
	_mm_storel_pd(&rL, radiusV);
	_mm_storeh_pd(&rH, radiusV);
	
	const __m128d thermV = _mm_set_pd((rH < heatingRadius) ? heatVal1 : heatVal2, // _mm_set_pd() is backwards
									  (rL < heatingRadius) ? heatVal1 : heatVal2)*randV;
	
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	
	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(cos(thetaH), cos(thetaL)));
}

void ThermalForceLocalized::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the localized thermal force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ThermalForce bit:
		forceFlags |= ThermalForceLocalizedFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
	{
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue1"), heatVal1, 6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue2"), heatVal2, 6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingRadius"), heatingRadius, 6, const_cast<char *> ("[m] (ThermalForceLocalized)"), error);
	}
}

void ThermalForceLocalized::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
	{
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue1"), &heatVal1, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue2"), &heatVal2, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingRadius"), &heatingRadius, NULL, error);
	}
}
