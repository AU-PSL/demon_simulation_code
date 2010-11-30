/*===- RectConfinementForce.cpp - libSimulation -===============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "RectConfinementForce.h"
	
RectConfinementForce::RectConfinementForce(Cloud *myCloud, double confineConstX, double confineConstY)
: Force(myCloud), confineX(-confineConstX), confineY(-confineConstY) {}

void RectConfinementForce::force1(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]), _mm_load_pd(&cloud->y[currentParticle]));
	}
}

void RectConfinementForce::force2(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l1[currentParticle])/v2, 
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n1[currentParticle])/v2);
	}
}

void RectConfinementForce::force3(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l2[currentParticle])/v2, 
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n2[currentParticle])/v2);
	}
}

void RectConfinementForce::force4(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l3[currentParticle]),
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n3[currentParticle]));
	}
}

inline void RectConfinementForce::force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY)
{
	const __m128d cVx = _mm_set1_pd(confineX);
	const __m128d cVy = _mm_set1_pd(confineY);
	double *pFx = &cloud->forceX[currentParticle];
	double *pFy = &cloud->forceY[currentParticle];

	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_mul_pd(cVx, currentPositionX));
	_mm_store_pd(pFy, _mm_load_pd(pFy) + _mm_mul_pd(cVy, currentPositionY));
}

void RectConfinementForce::writeForce(fitsfile *file, int *error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	//add flag indicating that the rectangular confinement force is used:
	if(!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add RectConfinementForce bit:
		forceFlags |= RectConfinementForceFlag;		//compound bitwise OR

		if(*error == 202 || *error == 204)	//keyword does not exist yet
			*error = 0;			//clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConstX"), confineX, 6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("confineConstY"), confineY, 6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
	}
}

void RectConfinementForce::readForce(fitsfile *file, int *error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if(!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConstX"), &confineX, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("confineConstY"), &confineY, NULL, error);
	}
}
