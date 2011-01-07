/*===- DrivingForce.cpp - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "DrivingForce.h"
#include <cmath>

const double DrivingForce::waveNum = 2.0*M_PI/0.002; //wavelength = 2mm
const double DrivingForce::angFreq = 2.0*M_PI*10.0; //10Hz

DrivingForce::DrivingForce(Cloud * const myCloud, const double drivingConst, const double amp, const double drivingShift)
: Force(myCloud), amplitude(amp), driveConst(-drivingConst), shift(drivingShift) {}

//1D:
void DrivingForce::force1_1D(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx1_pd(currentParticle));
}

void DrivingForce::force2_1D(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx2_pd(currentParticle));
}

void DrivingForce::force3_1D(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx3_pd(currentParticle));
}

void DrivingForce::force4_1D(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, vtime, cloud->getx4_pd(currentParticle));
}

//2D:
void DrivingForce::force1_2D(const double currentTime)
{
	force1_1D(currentTime); //DrivingForce affects only x-dimension
}

void DrivingForce::force2_2D(const double currentTime)
{
	force2_1D(currentTime);
}

void DrivingForce::force3_2D(const double currentTime)
{
	force3_1D(currentTime);
}

void DrivingForce::force4_2D(const double currentTime)
{
	force4_1D(currentTime);
}

//3D:
void DrivingForce::force1_3D(const double currentTime)
{
	force1_1D(currentTime); //DrivingForce affects only x-dimension
}

void DrivingForce::force2_3D(const double currentTime)
{
	force2_1D(currentTime);
}

void DrivingForce::force3_3D(const double currentTime)
{
	force3_1D(currentTime);
}

void DrivingForce::force4_3D(const double currentTime)
{
	force4_1D(currentTime);
}

//DrivingForce only acts in a single dimension, and privite force function is not inherited
// from Force.h, so no need to overload generic force function.
inline void DrivingForce::force(const unsigned int currentParticle, const __m128d currentTime, const __m128d currentPositionX)
{
	//N.B. F = A*sin(k*x - w*t)*exp(-(x + x0)^2/B) is the equation used in the paper, which differs from that below.
	const __m128d distV = currentPositionX - _mm_set1_pd(shift);
	const __m128d sinArg = _mm_set1_pd(waveNum)*currentPositionX - _mm_set1_pd(angFreq)*currentTime;
	const __m128d expArg = _mm_set1_pd(-1.0)*distV*distV/_mm_set1_pd(driveConst);

	//no SIMD trig instructions; break vectors and perform separately:
	double sinArgL, sinArgH, expArgL, expArgH;
	_mm_storel_pd(&sinArgL, sinArg);
	_mm_storeh_pd(&sinArgH, sinArg);
	_mm_storel_pd(&expArgL, expArg);
	_mm_storeh_pd(&expArgH, expArg);

	double * const pFx = cloud->forceX + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_set1_pd(amplitude)*_mm_set_pd(sin(sinArgH), sin(sinArgL))*_mm_set_pd(exp(expArgH), exp(expArgL))); // _mm_set_pd() is backwards
}

void DrivingForce::writeForce(fitsfile * const file, int * const error, const int dimension) const
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the driving force is used:
	if(!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add DrivingForce bit:
		forceFlags |= DrivingForceFlag; //compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;             //clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("drivingAmplitude"), amplitude, 6, const_cast<char *> ("[N] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingConst"), driveConst, 6, const_cast<char *> ("[m^2] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingShift"), shift, 6, const_cast<char *> ("[m] (DrivingForce)"), error);
	}
}

void DrivingForce::readForce(fitsfile * const file, int * const error, const int dimension)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if(!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("drivingAmplitude"), &amplitude, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingConst"), &driveConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingShift"), &shift, NULL, error);
	}
}
