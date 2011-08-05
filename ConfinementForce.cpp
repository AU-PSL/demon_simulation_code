/*===- ConfinementForce.cpp - libSimulation -===================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ConfinementForce.h"

//Constructors:
ConfinementForce1D::ConfinementForce1D(Cloud * const myCloud, double confineConst, double plasmaPotential) : Force(myCloud), confine(confineConst), potentialOffset(plasmaPotential) {}
ConfinementForce2D::ConfinementForce2D(Cloud * const myCloud, double confineConst, double plasmaPotential) : ConfinementForce1D(myCloud, confineConst, plasmaPotential) {}
ConfinementForce3D::ConfinementForce3D(Cloud * const myCloud, double confineConst, double plasmaPotential) : ConfinementForce2D(myCloud, confineConst, plasmaPotential) {}

//1D:
void ConfinementForce1D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void ConfinementForce1D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void ConfinementForce1D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void ConfinementForce1D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

//2D:
void ConfinementForce2D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void ConfinementForce2D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void ConfinementForce2D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void ConfinementForce2D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

//3D:
void ConfinementForce3D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), cloud->getz1_pd(currentParticle), cloud->getq1_pd(currentParticle));
}

void ConfinementForce3D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), cloud->getz2_pd(currentParticle), cloud->getq2_pd(currentParticle));
}

void ConfinementForce3D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), cloud->getz3_pd(currentParticle), cloud->getq3_pd(currentParticle));
}

void ConfinementForce3D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), cloud->getz4_pd(currentParticle), cloud->getq4_pd(currentParticle));
}

//force methods:
inline void ConfinementForce1D::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d charge)
{
	const __m128d cV = _mm_set1_pd(confine);
	double * const pFx = cloud->forceX + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + cV*charge*currentPositionX);

	//calculate contribution to electric potential:
	double * pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + _mm_set1_pd(-0.5)*cV*currentPositionX*currentPositionX + _mm_set1_pd(potentialOffset));
	
}

inline void ConfinementForce2D::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d charge)
{
	const __m128d cV = _mm_set1_pd(confine);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + cV*charge*currentPositionX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + cV*charge*currentPositionY);

	const __m128d rr = currentPositionX*currentPositionX + currentPositionY*currentPositionY;
	double * pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + _mm_set1_pd(-0.5)*cV*rr + _mm_set1_pd(potentialOffset));
}

inline void ConfinementForce3D::force(const cloud_index currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d currentPositionZ, const __m128d charge)
{
	const __m128d cV = _mm_set1_pd(confine);
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + cV*charge*currentPositionX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + cV*charge*currentPositionY);
	_mm_store_pd(pFz, _mm_load_pd(pFz) + cV*charge*currentPositionZ);

	const __m128d rr = currentPositionX*currentPositionX + currentPositionY*currentPositionY + currentPositionZ*currentPositionZ;
	double * pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + _mm_set1_pd(-0.5)*cV*rr + _mm_set1_pd(potentialOffset));
}

//write force:
void ConfinementForce1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	//add flag indicating that the confinement force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add ConfinementForce bit:
		forceFlags |= ConfinementForceFlag; //compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;                 //clear above error

		//add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("confineConst"), confine, 6, const_cast<char *> ("[N/m] (ConfinementForce)"), error);

	//write background plasma potential offset:
	if (!*error)
		fits_write_key_dbl(file, const_cast<char *> ("plasmaPotential"), potentialOffset, 6, const_cast<char *> ("[V] (background plasma potential offset)"), error);
}

//read force:
void ConfinementForce1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("confineConst"), &confine, NULL, error);
}
