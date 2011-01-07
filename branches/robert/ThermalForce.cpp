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

ThermalForce::ThermalForce(Cloud * const myCloud, const double redFactor) 
: Force(myCloud), mt(time(NULL)), heatVal(redFactor) {}

//1D:
void ThermalForce::force1_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle);
}

void ThermalForce::force2_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle);
}

void ThermalForce::force3_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle);
}

void ThermalForce::force4_1D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force1D(currentParticle);
}

//2D:
void ThermalForce::force1_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle);
}

void ThermalForce::force2_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle);
}

void ThermalForce::force3_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle);
}

void ThermalForce::force4_2D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force2D(currentParticle);
}

//3D:
void ThermalForce::force1_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle);
}

void ThermalForce::force2_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle);
}

void ThermalForce::force3_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle);
}

void ThermalForce::force4_3D(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force3D(currentParticle);
}

inline void ThermalForce::force1D(const unsigned int currentParticle)
{	
	//MT returns random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*_mm_set_pd(mt(), mt());
	const double twoL = mt()*2.0; //random number in (0,2)
	const double twoH = mt()*2.0;
	int directionL;
	int directionH;

	if(twoL < 1)
		directionL = -1; //left
	else if(twoL > 1)
		directionL = 1;  //right
	else if(twoL == 1)       //unlikely, but possible
		directionL = 0;  //no kick

	if(twoH < 1)
		directionH = -1;
	else if(twoH > 1)
		directionH = 1;
	else if(twoH == 1)
		directionH = 0;

	double * const pFx = cloud->forceX + currentParticle;
	
	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(directionH, directionL)); // _mm_set_pd() is backwards
}

inline void ThermalForce::force2D(const unsigned int currentParticle)
{	
	//MT returns random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*_mm_set_pd(mt(), mt());
	const double phiL = mt()*2.0*M_PI;	//azimuthal angle phi
	const double phiH = mt()*2.0*M_PI;
	
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	
	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(cos(phiH), cos(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(phiH), sin(phiL)));
}

inline void ThermalForce::force3D(const unsigned int currentParticle)
{	
	//MT returns random number in (0,1)
	const __m128d thermV = _mm_set1_pd(heatVal)*_mm_set_pd(mt(), mt());
	const double phiL = mt()*2.0*M_PI;	//azimuthal angle phi
	const double phiH = mt()*2.0*M_PI;
	const double thetaL = mt()*M_PI;	//polar angle theta
	const double thetaH = mt()*M_PI;
	
	double * const pFx = cloud->forceX + currentParticle;
	double * const pFy = cloud->forceY + currentParticle;
	double * const pFz = cloud->forceZ + currentParticle;
	
	_mm_store_pd(pFx, _mm_load_pd(pFx) + thermV*_mm_set_pd(sin(thetaH), sin(thetaL))*_mm_set_pd(cos(phiH), cos(phiL))); // _mm_set_pd() is backwards
	_mm_store_pd(pFy, _mm_load_pd(pFy) + thermV*_mm_set_pd(sin(thetaH), cos(thetaL))*_mm_set_pd(sin(phiH), sin(phiL)));
	_mm_store_pd(pFz, _mm_load_pd(pFz) + thermV*_mm_set_pd(cos(thetaH), cos(thetaL)));	
}

void ThermalForce::writeForce(fitsfile * const file, int * const error, const int dimension) const
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
		forceFlags |= ThermalForceFlag;		//compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;			//clear above error.

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

void ThermalForce::readForce(fitsfile * const file, int * const error, const int dimension)
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
