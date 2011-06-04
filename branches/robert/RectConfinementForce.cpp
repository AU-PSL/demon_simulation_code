/*===- RectConfinementForce.cpp - libSimulation -===============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "RectConfinementForce.h"
#include "ConfinementForce.h"     //used for 1D

//Constructors:
RectConfinementForce1D::RectConfinementForce1D(Cloud * const myCloud, double confineConstX)
: Force(myCloud), confineX(-confineConstX) {}
RectConfinementForce2D::RectConfinementForce2D(Cloud * const myCloud, double confineConstX, double confineConstY)
: RectConfinementForce1D(myCloud, confineConstX), confineY(-confineConstY) {}
RectConfinementForce3D::RectConfinementForce3D(Cloud * const myCloud, double confineConstX, double confineConstY, double confineConstZ)
: RectConfinementForce2D(myCloud, confineConstX, confineConstY), confineZ(-confineConstZ) {}

//1D:
void RectConfinementForce1D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle));
}

void RectConfinementForce1D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle));
}

void RectConfinementForce1D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle));
}

void RectConfinementForce1D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle));
}

//2D:
void RectConfinementForce2D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
}

void RectConfinementForce2D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
}

void RectConfinementForce2D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
}

void RectConfinementForce2D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
}

//3D:
void RectConfinementForce3D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), cloud->getz1_pd(currentParticle));
}

void RectConfinementForce3D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), cloud->getz2_pd(currentParticle));
}

void RectConfinementForce3D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), cloud->getz3_pd(currentParticle));
}

void RectConfinementForce3D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), cloud->getz4_pd(currentParticle));
}

//force methods:
inline void RectConfinementForce1D::force(const cloud_index currentParticle, const __m128d currentPositionX)
{
	double * const pFx = cloud->forceX + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_set1_pd(confineX)*currentPositionX);
}

inline void RectConfinementForce2D::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY)
{
	RectConfinementForce1D::force(currentParticle, currentPositionX);

	double * const pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFy, _mm_load_pd(pFy) + _mm_set1_pd(confineY)*currentPositionY);
}

inline void RectConfinementForce3D::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d currentPositionZ)
{
	RectConfinementForce2D::force(currentParticle, currentPositionX, currentPositionY);

	double * const pFz = cloud->forceZ + currentParticle;
	_mm_store_pd(pFz, _mm_load_pd(pFz) + _mm_set1_pd(confineZ)*currentPositionZ);
}

//writeForce:
void RectConfinementForce1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	//add flag indicating that the rectangular confinement force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add RectConfinementForce bit:
		forceFlags |= RectConfinementForceFlag; //compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;                     //clear above error.

		//add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConstX"), confineX, 6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
}

void RectConfinementForce2D::writeForce(fitsfile * const  file, int * const error) const
{
	RectConfinementForce1D::writeForce(file, error);
	fits_write_key_dbl(file, const_cast<char *> ("confineConstY"), confineY, 6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
}

void RectConfinementForce3D::writeForce(fitsfile * const  file, int * const error) const
{
	RectConfinementForce2D::writeForce(file, error);
	fits_write_key_dbl(file, const_cast<char *> ("confineConstZ"), confineZ, 6, const_cast<char *> ("[N/m] (RectConfinementForce)"), error);
}

//readForce:
void RectConfinementForce1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConstX"), &confineX, NULL, error);
}

void RectConfinementForce2D::readForce(fitsfile * const file, int * const error)
{
	RectConfinementForce1D::readForce(file, error);
	fits_read_key_dbl(file, const_cast<char *> ("confineConstY"), &confineY, NULL, error);
}

void RectConfinementForce3D::readForce(fitsfile * const file, int * const error)
{
	RectConfinementForce2D::readForce(file, error);
	fits_read_key_dbl(file, const_cast<char *> ("confineConstZ"), &confineZ, NULL, error);
}
