/*===- RotationalForce.cpp - libSimulation -====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "RotationalForce.h"
#include <cmath>

using namespace std;

RotationalForce::RotationalForce(Cloud *myCloud, double rmin, double rmax, double rotConst)
: Force(myCloud), rotationalConst(rotConst), innerRad(rmin), outerRad(rmax) {}

void RotationalForce::force1(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]), _mm_load_pd(&cloud->y[currentParticle]));
	}
}

void RotationalForce::force2(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2) 
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l1[currentParticle])/v2, 
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n1[currentParticle])/v2);
	}
}

void RotationalForce::force3(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l2[currentParticle])/v2, 
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n2[currentParticle])/v2);
	}
}

void RotationalForce::force4(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n; currentParticle < numParticles; currentParticle += 2)
	{
		force(currentParticle, _mm_load_pd(&cloud->x[currentParticle]) + _mm_load_pd(&cloud->l3[currentParticle]),
			_mm_load_pd(&cloud->y[currentParticle]) + _mm_load_pd(&cloud->n3[currentParticle]));
	}
}

inline void RotationalForce::force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY)
{
	const __m128d cIR = _mm_set1_pd(innerRad);
	const __m128d cOR = _mm_set1_pd(outerRad);
	const __m128d cR = _mm_sqrt_pd(currentPositionX*currentPositionX + currentPositionY*currentPositionY);

	// cR > cIR && cR < cOR
	const __m128d compV = _mm_and_pd(_mm_cmpgt_pd(cR, cIR), _mm_cmplt_pd(cR, cOR));
	double compL, compH;
	_mm_storel_pd(&compL, compV);
	_mm_storeh_pd(&compH, compV);
	
	const bool nanL = isnan(compL);
	const bool nanH = isnan(compH);
	if (!nanL && !nanH) //niether in, early return
		return;
	
	__m128d cRotConst = _mm_set_pd(nanH ? rotationalConst : 0.0, // _mm_set_pd() is backwards.
								   nanL ? rotationalConst : 0.0);
	
	//force in theta direction:
	double *pFx = &cloud->forceX[currentParticle];
	double *pFy = &cloud->forceY[currentParticle];
	
	// Fx = -c*x/r;
	// Fy = c*y/r;
	_mm_store_pd(pFx, _mm_load_pd(pFx) - cRotConst*currentPositionY/cR);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + cRotConst*currentPositionX/cR);
}

void RotationalForce::writeForce(fitsfile *file, int *error)
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

		//add RotationalForce bit:
		forceFlags |= RotationalForceFlag;		//compound bitwise OR

		if(*error == 202 && *error == 204)	//keyword does not exist yet
			*error = 0;			//clear above error.

		//add or update keyword:
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
	{
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("rotationalConst"), rotationalConst, 6, const_cast<char *> ("[N] (RotationalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("innerRadius"), innerRad, 6, const_cast<char *> ("[m] (RotationalForce)"), error);
		fits_write_key_dbl(file, const_cast<char *> ("outerRadius"), outerRad, 6, const_cast<char *> ("[m] (RotationalForce)"), error);
	}
}

void RotationalForce::readForce(fitsfile *file, int *error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if(!*error)
	{
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("rotationalConst"), &rotationalConst, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("innerRadius"), &innerRad, NULL, error);
		fits_read_key_dbl(file, const_cast<char *> ("outerRadius"), &outerRad, NULL, error);
	}
}
