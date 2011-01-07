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

//1D:
void DragForce::force1_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle, cloud->getVx1_pd(currentParticle));
}

void DragForce::force2_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle, cloud->getVx2_pd(currentParticle));
}

void DragForce::force3_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle, cloud->getVx3_pd(currentParticle));
}

void DragForce::force4_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle, cloud->getVx4_pd(currentParticle));
}

//2D:
void DragForce::force1_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle));
}

void DragForce::force2_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle));
}

void DragForce::force3_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle));
}

void DragForce::force4_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle));
}

//3D:
void DragForce::force1_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle), cloud->getVz1_pd(currentParticle));
}

void DragForce::force2_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle), cloud->getVz2_pd(currentParticle));
}

void DragForce::force3_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle), cloud->getVz3_pd(currentParticle));
}

void DragForce::force4_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle), cloud->getVz4_pd(currentParticle));
}

inline void DragForce::force1D(const unsigned int currentParticle, const __m128d currentVelocityX)
{
	const __m128d drag = _mm_set1_pd(dragConst)*_mm_load_pd(&cloud->mass[currentParticle]);
	double * const pFx = cloud->forceX + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + drag*currentVelocityX);
}

inline void DragForce::force2D(const unsigned int currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY)
{
	const __m128d drag = _mm_set1_pd(dragConst)*_mm_load_pd(&cloud->mass[currentParticle]);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + drag*currentVelocityX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + drag*currentVelocityY);
}

inline void DragForce::force3D(const unsigned int currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentVelocityZ)
{
	const __m128d drag = _mm_set1_pd(dragConst)*_mm_load_pd(&cloud->mass[currentParticle]);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + drag*currentVelocityX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + drag*currentVelocityY);
	_mm_store_pd(pFz, _mm_load_pd(pFz) + drag*currentVelocityZ);
}

void DragForce::writeForce(fitsfile * const file, int * const error, const int dimension) const
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
		forceFlags |= DragForceFlag; //compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;          //clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}
	
	if(!*error)
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("dragConst"), dragConst, 6, const_cast<char *> ("[s^-1] (DragForce)"), error);
}

void DragForce::readForce(fitsfile * const file, int * const error, const int dimension)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if(!*error)
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("dragConst"), &dragConst, NULL, error);
}
