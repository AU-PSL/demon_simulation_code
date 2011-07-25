/*===- MagneticForce.cpp - libSimulation -======================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "MagneticForce.h"

MagneticForce1D::MagneticForce1D(Cloud * const myCloud) : Force(myCloud) {}
MagneticForce2D::MagneticForce2D(Cloud * const myCloud, const double Bz) : MagneticForce1D(myCloud), Bz(Bz) {}
MagneticForce3D::MagneticForce3D(Cloud * const myCloud, const double Bz, const double By, const double Bx) : MagneticForce2D(myCloud, Bz), By(By), Bx(Bx) {}

void MagneticForce1D::force1(const double currentTime)
{
	// No magnetic force if particles restrained to 1D
}

void MagneticForce1D::force2(const double currentTime) {}
void MagneticForce1D::force3(const double currentTime) {}
void MagneticForce1D::force4(const double currentTime) {}

void MagneticForce2D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void MagneticForce2D::force2(const double currentTime)
{	
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void MagneticForce2D::force3(const double currentTime)
{	
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void MagneticForce2D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

void MagneticForce3D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx1_pd(currentParticle), cloud->getVy1_pd(currentParticle), cloud->getVz1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void MagneticForce3D::force2(const double currentTime)
{	
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx2_pd(currentParticle), cloud->getVy2_pd(currentParticle), cloud->getVz2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void MagneticForce3D::force3(const double currentTime)
{	
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx3_pd(currentParticle), cloud->getVy3_pd(currentParticle), cloud->getVz3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void MagneticForce3D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getVx4_pd(currentParticle), cloud->getVy4_pd(currentParticle), cloud->getVz4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

inline void MagneticForce2D::force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentCharge)
{
	const __m128d qBz = currentCharge*_mm_set1_pd(Bz);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + qBz*currentVelocityY);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - qBz*currentVelocityX);
}

inline void MagneticForce3D::force(const cloud_index currentParticle, const __m128d currentVelocityX, const __m128d currentVelocityY, const __m128d currentVelocityZ, const __m128d currentCharge)
{
	const __m128d qBx = currentCharge*_mm_set1_pd(Bx);
	const __m128d qBy = currentCharge*_mm_set1_pd(By);
	const __m128d qBz = currentCharge*_mm_set1_pd(Bz);

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx)  + qBz*currentVelocityY - qBy*currentVelocityZ);
	_mm_store_pd(pFy, _mm_load_pd(pFy)  + qBx*currentVelocityZ - qBz*currentVelocityZ);
	_mm_store_pd(pFz, _mm_load_pd(pFz)  + qBy*currentVelocityX - qBx*currentVelocityY);
}

void MagneticForce1D::writeForce(fitsfile * const file, int * const error) const {}

void MagneticForce2D::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the magnetic force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add MagneticForce bit:
		forceFlags |= MagneticForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}
	
	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("magnetic field (z)"), Bz, 6, const_cast<char *> ("[T] (MagneticForce)"), error);
}

void MagneticForce3D::writeForce(fitsfile * const file, int * const error) const
{
	MagneticForce2D::writeForce(file, error);
	if (!*error)
		fits_write_key_dbl(file, const_cast<char *> ("magnetic field (y)"), By, 6, const_cast<char *> ("[T] (MagneticForce)"), error);

	if (!*error)
		fits_write_key_dbl(file, const_cast<char *> ("magnetic field (x)"), Bx, 6, const_cast<char *> ("[T] (MagneticForce)"), error);
}

void MagneticForce1D::readForce(fitsfile * const file, int * const error) {}

void MagneticForce2D::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("magnetic field (z)"), &Bz, NULL, error);
}

void MagneticForce3D::readForce(fitsfile * const file, int * const error)
{
	MagneticForce2D::readForce(file, error);

	if (!*error)
		fits_read_key_dbl(file, const_cast<char *> ("magnetic field (y)"), &By, NULL, error);

	if (!*error)
		fits_read_key_dbl(file, const_cast<char *> ("magnetic field (x)"), &Bx, NULL, error);
}
