/*===- ThermalForceLocalized.cpp - libSimulation -==============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#include "ThermalForceLocalized.h"
#include <cmath>
#include <ctime>

ThermalForceLocalized1D::ThermalForceLocalized1D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius) 
: Force(myCloud), mt(time(NULL)), heatingRadius(specifiedRadius), heatVal1(thermRed1), heatVal2(thermRed2) {}
ThermalForceLocalized2D::ThermalForceLocalized2D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius) 
: ThermalForceLocalized1D(myCloud, thermRed1, thermRed2, specifiedRadius) {}
ThermalForceLocalized3D::ThermalForceLocalized3D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius) 
: ThermalForceLocalized2D(myCloud, thermRed1, thermRed2, specifiedRadius) {}

//1D:
void ThermalForceLocalized1D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle));
}

void ThermalForceLocalized1D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle));
}

void ThermalForceLocalized1D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx3_pd(currentParticle));
}

void ThermalForceLocalized1D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle));
}

//2D:
void ThermalForceLocalized2D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle));
}

void ThermalForceLocalized2D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle));
}

void ThermalForceLocalized2D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle));
}

void ThermalForceLocalized2D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle));
}

//3D:
void ThermalForceLocalized3D::force1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx1_pd(currentParticle), cloud->gety1_pd(currentParticle), cloud->getz1_pd(currentParticle));
}

void ThermalForceLocalized3D::force2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx2_pd(currentParticle), cloud->gety2_pd(currentParticle), cloud->getz2_pd(currentParticle));
}

void ThermalForceLocalized3D::force3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx3_pd(currentParticle), cloud->gety3_pd(currentParticle), cloud->getz3_pd(currentParticle));
}

void ThermalForceLocalized3D::force4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, cloud->getx4_pd(currentParticle), cloud->gety4_pd(currentParticle), cloud->getz4_pd(currentParticle));
}

//Force methods:
inline void ThermalForceLocalized1D::force(const cloud_index currentParticle, const __m128d displacementX)
{
	const __m128d radiusV = _mm_sqrt_pd(displacementX*displacementX); //absolute value

	//MT returns random number in (0,1)
	const double directionL = mt()*2.0 - 1.0; //random kick
	const double directionH = mt()*2.0 - 1.0; //random number in (-1, 1)

	double rL, rH;
	_mm_storel_pd(&rL, radiusV);
	_mm_storeh_pd(&rH, radiusV);

	const __m128d thermV = _mm_set_pd((rH < heatingRadius) ? heatVal1 : heatVal2, // _mm_set_pd() is backwards
		(rL < heatingRadius) ? heatVal1 : heatVal2);

	double * const pFx = cloud->forceX + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(directionH, directionL)); // _mm_set_pd() is backwards
}

inline void ThermalForceLocalized2D::force(const cloud_index currentParticle, const __m128d displacementX, const __m128d displacementY)
{
	const __m128d radiusV = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);

	const double phiL = mt()*2.0*M_PI;
	const double phiH = mt()*2.0*M_PI;

	double rL, rH;
	_mm_storel_pd(&rL, radiusV);
	_mm_storeh_pd(&rH, radiusV);

	const __m128d thermV = _mm_set_pd((rH < heatingRadius) ? heatVal1 : heatVal2, // _mm_set_pd() is backwards
		(rL < heatingRadius) ? heatVal1 : heatVal2)*_mm_set_pd(mt(), mt());

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(phiH), sin(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(phiH), sin(phiL)));
}

inline void ThermalForceLocalized3D::force(const cloud_index currentParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ)
{
	const __m128d radiusV = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY + displacementZ*displacementZ);

	const double phiL = mt()*2.0*M_PI;
	const double phiH = mt()*2.0*M_PI;
	const double thetaL = mt()*M_PI;
	const double thetaH = mt()*M_PI;

	double rL, rH;
	_mm_storel_pd(&rL, radiusV);
	_mm_storeh_pd(&rH, radiusV);

	const __m128d thermV = _mm_set_pd((rH < heatingRadius) ? heatVal1 : heatVal2, // _mm_set_pd() is backwards
		(rL < heatingRadius) ? heatVal1 : heatVal2)*_mm_set_pd(mt(), mt());

	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;

	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))*_mm_set_pd(cos(phiH), cos(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))*_mm_set_pd(sin(phiH), sin(phiL)));
	_mm_store_pd(pFz, _mm_load_pd(pFz) + thermV*_mm_set_pd(cos(thetaH), sin(thetaL)));
}

//writeForce:
void ThermalForceLocalized1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the localized thermal force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add ThermalForce bit:
		forceFlags |= ThermalForceLocalizedFlag; //compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;                     //clear above error.

		//add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue1"), heatVal1, 6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingValue2"), heatVal2, 6, const_cast<char *> ("[N] (ThermalForceLocalized)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("heatingRadius"), heatingRadius, 6, const_cast<char *> ("[m] (ThermalForceLocalized)"), error);
	}
}

//readForce:
void ThermalForceLocalized1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue1"), &heatVal1, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingValue2"), &heatVal2, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("heatingRadius"), &heatingRadius, NULL, error);
	}
}
