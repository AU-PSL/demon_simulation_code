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

const double ShieldedCoulombForce::coulomb = 1.0/(4.0*M_PI*Cloud::epsilon0);

ShieldedCoulombForce::ShieldedCoulombForce(Cloud * const myCloud, const double shieldingConstant)
: Force(myCloud), shielding(shieldingConstant) {}

void ShieldedCoulombForce::force1(const double currentTime)
{
    (void)currentTime;
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx1_pd(currentParticle);
		const __m128d vy1 = cloud->gety1_pd(currentParticle);
		const __m128d vq1 = cloud->getq1_pd(currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);

		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			force(currentParticle, i, vq1, cloud->getq1_pd(i), vx1 - cloud->getx1_pd(i), vy1 - cloud->gety1_pd(i));
			forcer(currentParticle, i, vq1, cloud->getq1r_pd(i), vx1 - cloud->getx1r_pd(i), vy1 - cloud->gety1r_pd(i));
		}
	}
}

void ShieldedCoulombForce::force2(const double currentTime)
{
    (void)currentTime;
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx2_pd(currentParticle);
		const __m128d vy1 = cloud->gety2_pd(currentParticle);
		const __m128d vq1 = cloud->getq2_pd(currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			force(currentParticle, i, vq1, cloud->getq2_pd(i), vx1 - cloud->getx2_pd(i), vy1 - cloud->gety2_pd(i));
			forcer(currentParticle, i, vq1, cloud->getq2r_pd(i), vx1 - cloud->getx2r_pd(i), vy1 - cloud->gety2r_pd(i));
		}
	}
}

void ShieldedCoulombForce::force3(const double currentTime)
{
    (void)currentTime;
    for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx3_pd(currentParticle);
		const __m128d vy1 = cloud->gety3_pd(currentParticle);
		const __m128d vq1 = cloud->getq3_pd(currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			force(currentParticle, i, vq1, cloud->getq3_pd(i), vx1 - cloud->getx3_pd(i), vy1 - cloud->gety3_pd(i));
			forcer(currentParticle, i, vq1, cloud->getq3r_pd(i), vx1 - cloud->getx3r_pd(i), vy1 - cloud->gety3r_pd(i));
		}
	}
}

void ShieldedCoulombForce::force4(const double currentTime)
{
    (void)currentTime;
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx4_pd(currentParticle);
		const __m128d vy1 = cloud->gety4_pd(currentParticle);
		const __m128d vq1 = cloud->getq4_pd(currentParticle);
		double x1, x2, y1, y2, q1, q2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		_mm_storel_pd(&q1, vq1);
		_mm_storeh_pd(&q2, vq1);

		force(currentParticle, currentParticle + 1, q1, q2, x1 - x2, y1 - y2);
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			force(currentParticle, i, vq1, cloud->getq4_pd(i), vx1 - cloud->getx4_pd(i), vy1 - cloud->gety4_pd(i));
			forcer(currentParticle, i, vq1, cloud->getq4r_pd(i), vx1 - cloud->getx4r_pd(i), vy1 - cloud->gety4r_pd(i));
		}
	}
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const double currentCharge, const double iCharge,
                                        const double displacementX, const double displacementY)
{
	// Calculate displacement between particles.
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement*shielding;

	if (valExp < 10.0) // restrict to 10*(ion debye length)
	{
		// calculate phi
		const double coefficient = coulomb/(displacement*exp(valExp));
		cloud->phi[currentParticle] += coefficient*iCharge;
		cloud->phi[iParticle] += coefficient*currentCharge;

		// calculate force
		const double forceC = currentCharge*iCharge*coefficient*(1.0 + valExp)/(displacement*displacement);
		cloud->forceX[currentParticle] += forceC*displacementX;
		cloud->forceY[currentParticle] += forceC*displacementY;

		// equal and opposite force:
		cloud->forceX[iParticle] -= forceC*displacementX;
		cloud->forceY[iParticle] -= forceC*displacementY;
	}
}

inline void ShieldedCoulombForce::force(const cloud_index currentParticle, const cloud_index iParticle, 
                                        const __m128d currentCharge, const __m128d iCharge,
                                        const __m128d displacementX, const __m128d displacementY)
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
	
	__m128d expv = _mm_set_pd(boolH ? exp(-valExpH) : 0.0, // _mm_set_pd is backwards
							  boolL ? exp(-valExpL) : 0.0);

	// calculate phi
	const __m128d coefficient = _mm_set1_pd(coulomb)/displacement*expv;
	double *pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + coefficient*iCharge);
	pPhi = cloud->phi + iParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + coefficient*currentCharge);
	
	// calculate force
	const __m128d forceC = currentCharge*iCharge*coefficient*(_mm_set1_pd(1.0) + valExp)/(displacement*displacement);
	
	const __m128d forcevX = forceC*displacementX;
	const __m128d forcevY = forceC*displacementY;

	double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);

	// equal and opposite force:
	pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) - forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) - forcevY);
}

inline void ShieldedCoulombForce::forcer(const cloud_index currentParticle, const cloud_index iParticle, 
                                         const __m128d currentCharge, const __m128d iCharge,
                                         const __m128d displacementX, const __m128d displacementY)
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
	
	__m128d expv = _mm_set_pd(boolH ? exp(-valExpH) : 0.0, // _mm_set_pd is backwards
							  boolL ? exp(-valExpL) : 0.0);
	
	// calculate phi
	const __m128d coefficient = _mm_set1_pd(coulomb)/displacement*expv;
	double *pPhi = cloud->phi + currentParticle;
	_mm_store_pd(pPhi, _mm_load_pd(pPhi) + coefficient*iCharge);
	pPhi = cloud->phi + iParticle;
	_mm_storer_pd(pPhi, _mm_loadr_pd(pPhi) + coefficient*currentCharge);
    
	// calculate force
	const __m128d forceC = currentCharge*iCharge*coefficient*(_mm_set1_pd(1.0) + valExp)/(displacement*displacement);
	
	const __m128d forcevX = forceC*displacementX;
	const __m128d forcevY = forceC*displacementY;
	
	double *pFx = cloud->forceX + currentParticle;
	double *pFy = cloud->forceY + currentParticle;
	_mm_store_pd(pFx, _mm_load_pd(pFx) + forcevX);
	_mm_store_pd(pFy, _mm_load_pd(pFy) + forcevY);
	
	// equal and opposite force:
	pFx = cloud->forceX + iParticle;
	pFy = cloud->forceY + iParticle;
	_mm_storer_pd(pFx, _mm_loadr_pd(pFx) - forcevX);
	_mm_storer_pd(pFy, _mm_loadr_pd(pFy) - forcevY);
}

void ShieldedCoulombForce::writeForce(fitsfile * const file, int * const error) const
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	// add flag indicating that the shielded Coulomb force is used:
	if (!*error) 
	{
		long forceFlags = 0;
		fits_read_key_lng(file, const_cast<char *> ("FORCES"), &forceFlags, NULL, error);

		// add ShieldedCoulombForce bit:
		forceFlags |= ShieldedCoulombForceFlag; // compound bitwise OR

		if (*error == KEY_NO_EXIST || *error == VALUE_UNDEFINED)
			*error = 0; // clear above error.

		// add or update keyword.
		if (!*error) 
			fits_update_key(file, TLONG, const_cast<char *> ("FORCES"), &forceFlags, const_cast<char *> ("Force configuration."), error);
	}

	if (!*error)
		// file, key name, value, precision (scientific format), comment
		fits_write_key_dbl(file, const_cast<char *> ("shieldingConstant"), shielding, 6, const_cast<char *> ("[m^-1] (ShieldedCoulombForce)"), error);
}

void ShieldedCoulombForce::readForce(fitsfile * const file, int * const error)
{
	// move to primary HDU:
	if (!*error)
		// file, # indicating primary HDU, HDU type, error
 		fits_movabs_hdu(file, 1, IMAGE_HDU, error);
	
	if (!*error)
		// file, key name, value, don't read comment, error
		fits_read_key_dbl(file, const_cast<char *> ("shieldingConstant"), &shielding, NULL, error);
}
