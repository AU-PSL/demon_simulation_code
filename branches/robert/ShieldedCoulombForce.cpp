/*===- ShieldedCoulombForce.cpp - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "ShieldedCoulombForce.h"
#include <cmath>

ShieldedCoulombForce::ShieldedCoulombForce(Cloud * const myCloud, const double shieldingConstant)
: Force(myCloud), shielding(shieldingConstant) {}

void ShieldedCoulombForce::force1(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx1_pd(currentParticle);
		const __m128d vy1 = cloud->gety1_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);

		force(currentParticle, currentParticle + 1, x1 - x2, y1 - y2);

		for(unsigned int i = currentParticle + 2; i < numParticles; i += 2)
		{
			const double * const px2 = cloud->x + i;	//increment memory location
			const double * const py2 = cloud->y + i;

			force(currentParticle, i, vx1 - _mm_load_pd(px2), vy1 - _mm_load_pd(py2));
			forcer(currentParticle, i, vx1 - _mm_loadr_pd(px2), vy1 - _mm_loadr_pd(py2));
		}
	}
}

void ShieldedCoulombForce::force2(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx2_pd(currentParticle);
		const __m128d vy1 = cloud->gety2_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);

		force(currentParticle, currentParticle + 1, x1 - x2, y1 - y2);
		for(unsigned int i = currentParticle + 2; i < numParticles; i += 2)
		{
			const double * const px2 = cloud->x + i;	//increment memory location
			const double * const py2 = cloud->y + i;
			const double * const pl = cloud->l1 + i;
			const double * const pn = cloud->n1 + i;
            
			force(currentParticle, i, vx1 - (_mm_load_pd(px2) + _mm_load_pd(pl)/v2), 
				vy1 - (_mm_load_pd(py2) + _mm_load_pd(pn)/v2));
            
			forcer(currentParticle, i, vx1 - (_mm_loadr_pd(px2) + _mm_loadr_pd(pl)/v2), 
				vy1 - (_mm_loadr_pd(py2) + _mm_loadr_pd(pn)/v2));
		}
	}
}

void ShieldedCoulombForce::force3(const double currentTime)
{
	const __m128d v2 = _mm_set1_pd(2.0);
	for (unsigned int currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx3_pd(currentParticle);
		const __m128d vy1 = cloud->gety3_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);

		force(currentParticle, currentParticle + 1, x1 - x2, y1 - y2);
		for(unsigned int i = currentParticle + 2; i < numParticles; i += 2)
		{
			const double * const px2 = cloud->x + i;	//increment memory location
			const double * const py2 = cloud->y + i;
			const double * const pl = cloud->l2 + i;
			const double * const pn = cloud->n2 + i;

			force(currentParticle, i, vx1 - (_mm_load_pd(px2) + _mm_load_pd(pl)/v2), 
				vy1 - (_mm_load_pd(py2) + _mm_load_pd(pn)/v2));
            
			forcer(currentParticle, i, vx1 - (_mm_loadr_pd(px2) + _mm_loadr_pd(pl)/v2), 
 				vy1 - (_mm_loadr_pd(py2) + _mm_loadr_pd(pn)/v2));
		}
	}
}

void ShieldedCoulombForce::force4(const double currentTime)
{
	for (unsigned int currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx4_pd(currentParticle);
		const __m128d vy1 = cloud->gety4_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);

		force(currentParticle, currentParticle + 1, x1 - x2, y1 - y2);
		for(unsigned int i = currentParticle + 2; i < numParticles; i += 2)
		{
			const double * const px2 = cloud->x + i;	//increment memory location
			const double * const py2 = cloud->y + i;
			const double * const pl = cloud->l3 + i;
			const double * const pn = cloud->n3 + i;
            
			force(currentParticle, i, vx1 - (_mm_load_pd(px2) + _mm_load_pd(pl)), 
				vy1 - (_mm_load_pd(py2) + _mm_load_pd(pn)));

			forcer(currentParticle, i, vx1 - (_mm_loadr_pd(px2) + _mm_loadr_pd(pl)), 
 				vy1 - (_mm_loadr_pd(py2) + _mm_loadr_pd(pn)));
		}
	}
}

inline void ShieldedCoulombForce::force(const unsigned int currentParticle, const unsigned int iParticle, const double displacementX, const double displacementY)
{
	// Calculate displacement between particles.
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement*shielding;

	if(valExp < 10.0) 	//restrict to 10*(ion debye length)
	{
		 //conclude force calculation:
		const double displacement3 = displacement*displacement*displacement;
		//set to charges multiplied by Coulomb's constant:
		const double exponential = (cloud->charge[currentParticle]*cloud->charge[iParticle])/(4.0*M_PI*8.85E-12)*(1.0 + valExp)/(displacement3*exp(valExp));
		cloud->forceX[currentParticle] += exponential*displacementX;
		cloud->forceY[currentParticle] += exponential*displacementY;

		//equal and opposite force:
		cloud->forceX[iParticle] -= exponential*displacementX;
		cloud->forceY[iParticle] -= exponential*displacementY;
	}
}

inline void ShieldedCoulombForce::force(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY)
{
	// Calculate displacement between particles.
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement*_mm_set_pd(shielding, shielding);
	
	double valExpL, valExpH;
	_mm_storel_pd(&valExpL, valExp);
	_mm_storeh_pd(&valExpH, valExp);
	
	const bool boolL = valExpL < 10.0;
	const bool boolH = valExpH < 10.0;
	if (!boolL && !boolH)
		return;

	double expL, expH;
	_mm_storel_pd(&expL, valExp);
	_mm_storeh_pd(&expH, valExp);
	
	__m128d expv = _mm_set_pd(boolH ? exp(-expH) : 0.0, //_mm_set_pd is backwards
							  boolL ? exp(-expL) : 0.0);

	//conclude force calculation:
	const __m128d displacement3 = displacement*displacement*displacement;
	//set to charges multiplied by Coulomb's constant:
	const double c = 4.0*M_PI*8.85E-12;
	const __m128d exponential = _mm_load_pd(&cloud->charge[currentParticle])*_mm_load_pd(&cloud->charge[iParticle])
		/_mm_set_pd(c, c)*(_mm_set1_pd(1.0) + valExp)/displacement3*expv;
	
	const __m128d forcevX = exponential*displacementX;
	const __m128d forcevY = exponential*displacementY;

	double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);

	//equal and opposite force:
	pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) - forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - forcevY);
}

inline void ShieldedCoulombForce::forcer(const unsigned int currentParticle, const unsigned int iParticle, const __m128d displacementX, const __m128d displacementY)
{
	// Calculate displacement between particles.
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement*_mm_set_pd(shielding, shielding);
	
	double valExpL, valExpH;
	_mm_storel_pd(&valExpL, valExp);
	_mm_storeh_pd(&valExpH, valExp);
	
	const bool boolL = valExpL < 10.0;
	const bool boolH = valExpH < 10.0;
	if (!boolL && !boolH)
		return;
	
	double expL, expH;
	_mm_storel_pd(&expL, valExp);
	_mm_storeh_pd(&expH, valExp);
	
	__m128d expv = _mm_set_pd(boolH ? exp(-expH) : 0.0, //_mm_set_pd is backwards
							  boolL ? exp(-expL) : 0.0);
    
	//conclude force calculation:
	const __m128d displacement3 = displacement*displacement*displacement;
	//set to charges multiplied by Coulomb's constant:
	const double c = 4.0*M_PI*8.85e-12;
	const __m128d exponential = _mm_load_pd(&cloud->charge[currentParticle])*_mm_loadr_pd(&cloud->charge[iParticle])
		/_mm_set_pd(c, c)*(_mm_set1_pd(1.0) + valExp)/displacement3*expv;

	const __m128d forcevX = exponential*displacementX;
	const __m128d forcevY = exponential*displacementY;

	double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);

	//equal and opposite force:
	pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle; 
	_mm_storer_pd(pFx, _mm_loadr_pd(pFx) - forcevX);
	_mm_storer_pd(pFy, _mm_loadr_pd(pFy) - forcevY);
}

void ShieldedCoulombForce::writeForce(fitsfile * const file, int * const error) const
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

		//add ShieldedCoulombForce bit:
		forceFlags |= ShieldedCoulombForceFlag;		//compound bitwise OR

		if(*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0;			//clear above error.

		//add or update keyword.
		if(!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if(!*error)
		//file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("shieldingConstant"), shielding, 6, const_cast<char *> ("[m^-1] (ShieldedCoulombForce)"), error);
}

void ShieldedCoulombForce::readForce(fitsfile * const file, int * const error)
{
	//move to primary HDU:
	if(!*error)
		//file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if(!*error)
		//file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("shieldingConstant"), &shielding, NULL, error);
}
