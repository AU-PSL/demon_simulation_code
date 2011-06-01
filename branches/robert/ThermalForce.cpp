/*===- ThermalForce.cpp - libSimulation -=======================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "ThermalForce.h"

//Constructors:
ThermalForce1D::ThermalForce1D(Cloud * const myCloud, const double redFactor) : Force(myCloud), mt(time(NULL)), heatVal(redFactor) {}
ThermalForce2D::ThermalForce2D(Cloud * const myCloud, const double redFactor) : ThermalForce1D(myCloud, redFactor) {}
ThermalForce3D::ThermalForce3D(Cloud * const myCloud, const double redFactor) : ThermalForce2D(myCloud, redFactor) {}

//1D:
void ThermalForce1D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce1D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce1D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce1D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

//2D:
void ThermalForce2D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce2D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce2D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce2D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

//3D:
void ThermalForce3D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce3D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce3D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

void ThermalForce3D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle);
}

//force methods:
inline void ThermalForce1D::force(const cloud_index currentParticle)
{
	//MT returns random number in (0,1)
	const double directionL = mt()*2.0 - 1.0; //random kick
	const double directionH = mt()*2.0 - 1.0; //random number in (-1, 1)

	double * const pFx = cloud->forceX + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_set_pd(directionH, directionL)); // _mm_set_pd() is backwards
}

inline void ThermalForce2D::force(const cloud_index currentParticle)
{
	//MT returns random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*_mm_set_pd(mt(), mt()); //random strength
	const double phiL = mt()*2.0*M_PI; //random direction
	const double phiH = mt()*2.0*M_PI;

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(cos(phiH), cos(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(phiH), sin(phiL)));
}

inline void ThermalForce3D::force(const cloud_index currentParticle)
{
	//MT returns random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*_mm_set_pd(mt(), mt()); //random strength
	const double phiL = mt()*2.0*M_PI; //random azimuthal angle phi
	const double phiH = mt()*2.0*M_PI;
	const double thetaL = mt()*M_PI;   //random polar angle theta
	const double thetaH = mt()*M_PI;

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))*_mm_set_pd(cos(phiH), cos(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(thetaH), cos(thetaL))*_mm_set_pd(sin(phiH), sin(phiL)));
	_mm_store_pd(pFz, _mm_load_pd(pFz) + thermV*_mm_set_pd(cos(thetaH), cos(thetaL)));
}

//writeForce:
void ThermalForce1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the thermal force is used:
	if(!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add ThermalForce bit:
		forceFlags |= ThermalForceFlag; //compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;             //clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue"), heatVal, 6, const_cast<char *> ("[N] (ThermalForce)"), error);
	}
}

void ThermalForce1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if(!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue"), &heatVal, NULL, error);
	}
}
