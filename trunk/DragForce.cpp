/*===- DragForce.cpp - libSimulation -==========================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "DragForce.h"

DragForce::DragForce(Cloud * const myCloud, const double gamma) 
: Force(myCloud), dragConst(-gamma) {}

void DragForce::force1(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, _mm_load_pd(&cloud->Vx[currentParticle]), _mm_load_pd(&cloud->Vy[currentParticle]));
}

void DragForce::force2(const double currentTime)
{	
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, _mm_load_pd(&cloud->Vx[currentParticle]) + _mm_load_pd(&cloud->k1[currentParticle])/v2, 
			_mm_load_pd(&cloud->Vy[currentParticle]) + _mm_load_pd(&cloud->m1[currentParticle])/v2);
}

void DragForce::force3(const double currentTime)
{	
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, _mm_load_pd(&cloud->Vx[currentParticle]) + _mm_load_pd(&cloud->k2[currentParticle])/v2, 
			_mm_load_pd(&cloud->Vy[currentParticle]) + _mm_load_pd(&cloud->m2[currentParticle])/v2);
}

void DragForce::force4(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, _mm_load_pd(&cloud->Vx[currentParticle]) + _mm_load_pd(&cloud->k3[currentParticle]), 
			_mm_load_pd(&cloud->Vy[currentParticle]) + _mm_load_pd(&cloud->m3[currentParticle]));
}

inline void DragForce::force(const unsigned int currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY)
{
	const __m128d drag = _mm_set1_pd(dragConst)*_mm_load_pd(&cloud->mass[currentParticle]);
	double * const pFx = &cloud->forceX[currentParticle];
	double * const pFy = &cloud->forceY[currentParticle];

	_mm_store_pd(pFx, _mm_load_pd(pFx) + drag*currentVelocityX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + drag*currentVelocityY);
}

void DragForce::writeForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	//add flag indicating that the drag force is used:
	if(!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add DragForce bit:
		forceFlags |= DragForceFlag;		//compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;			//clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}
	
	if(!*error)
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("dragConst"), dragConst, 6, const_cast<char *> ("[s^-1] (DragForce)"), error);
}

void DragForce::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if(!*error)
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("dragConst"), &dragConst, NULL, error);
}
