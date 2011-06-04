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

const double DrivingForce1D::waveNum = 2.0*M_PI/0.002; //wavelength = 2mm
const double DrivingForce1D::angFreq = 2.0*M_PI*10.0; //10Hz

//Constructors:
DrivingForce1D::DrivingForce1D(Cloud * const myCloud, const double drivingConst, const double amp, const double drivingShift)
: Force(myCloud), amplitude(amp), driveConst(-drivingConst), shift(drivingShift) {}
DrivingForce2D::DrivingForce2D(Cloud * const myCloud, const double drivingConst, const double amp, const double drivingShift)
: DrivingForce1D(myCloud, drivingConst, amp, drivingShift) {}
DrivingForce3D::DrivingForce3D(Cloud * const myCloud, const double drivingConst, const double amp, const double drivingShift)
: DrivingForce2D(myCloud, drivingConst, amp, drivingShift) {}

//1D:
void DrivingForce1D::force1(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx1_pd(currentParticle));
}

void DrivingForce1D::force2(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx2_pd(currentParticle));
}

void DrivingForce1D::force3(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
		force(currentParticle, vtime, cloud->getx3_pd(currentParticle));
}

void DrivingForce1D::force4(const double currentTime)
{
	const __m128d vtime = _mm_set1_pd(currentTime);
	for (cloud_index currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
		force(currentParticle, vtime, cloud->getx4_pd(currentParticle));
}

//force method:
inline void DrivingForce1D::force(const cloud_index currentParticle, const __m128d currentTime, const __m128d currentPositionX)
{
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
	_mm_store_pd(pFx, _mm_load_pd(pFx) + _mm_set1_pd(amplitude)*_mm_set_pd(sin(sinArgH),
		sin(sinArgL))*_mm_set_pd(exp(expArgH), exp(expArgL))); // _mm_set_pd() is backwards
}

//writeForce:
void DrivingForce1D::writeForce(fitsfile * const file, int * const error) const
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	//add flag indicating that the driving force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		//add DrivingForce bit:
		forceFlags |= DrivingForceFlag; //compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;             //clear above error.

		//add or update keyword:
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("drivingAmplitude"), amplitude, 6, const_cast<char *> ("[N] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingConst"), driveConst, 6, const_cast<char *> ("[m^2] (DrivingForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("drivingShift"), shift, 6, const_cast<char *> ("[m] (DrivingForce)"), error);
	}
}

//readForce:
void DrivingForce1D::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if (!*error)
		//file, # indicating primary HDU, HDU type, error
		fits_movabs_hdu(file, 1, IMAGE_HDU, error);

	if (!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("drivingAmplitude"), &amplitude, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingConst"), &driveConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("drivingShift"), &shift, NULL, error);
	}
}
