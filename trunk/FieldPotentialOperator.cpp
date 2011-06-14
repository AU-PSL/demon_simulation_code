/*===- FieldPotentialOperator.cpp - libSimulation -=============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#include "FieldPotentialOperator.h"
#include <cmath>

void FieldPotentialOperator::operation1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d x1v = cloud->getx1_pd(currentParticle);
		const __m128d y1v = cloud->gety1_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, x1v);
		_mm_storeh_pd(&x2, x1v);
		_mm_storel_pd(&y1, y1v);
		_mm_storeh_pd(&y2, y1v);

		const double dispX = x1 - x2;
		const double dispY = y1 - y2;
		const double disp = sqrt(dispX*dispX + dispY*dispY);
		__m128d dispXv = _mm_set_pd(dispX, dispX);
		__m128d dispYv = _mm_set_pd(dispY, dispY);
		__m128d dispV = _mm_set_pd(disp, disp);

		potential(currentParticle, disp, cloud->getq1r_pd(currentParticle));
		field(currentParticle, dispV, dispXv, dispYv, cloud->getphi1_pd(currentParticle));

		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			dispXv = x1v - cloud->getx1_pd(i);
			dispYv = y1v - cloud->gety1_pd(i);
			dispV = _mm_sqrt_pd(dispXv*dispXv + dispYv*dispYv);
			const __m128d dispXvr = x1v - cloud->getx1r_pd(i);
			const __m128d dispYvr = y1v - cloud->gety1r_pd(i);
			const __m128d dispVr = _mm_sqrt_pd(dispXvr*dispXvr + dispYvr*dispYvr);

			potential(currentParticle, i, dispV, cloud->getq1_pd(currentParticle), cloud->getq1_pd(i));
			potentialr(currentParticle, i, dispVr, cloud->getq1_pd(currentParticle), cloud->getq1r_pd(i));

			field(currentParticle, i, dispV, dispXv, dispYv, cloud->getphi1_pd(currentParticle), cloud->getphi1_pd(i));
			fieldr(currentParticle, i, dispVr, dispXvr, dispYvr, cloud->getphi1_pd(currentParticle), cloud->getphi1r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d x2v = cloud->getx2_pd(currentParticle);
		const __m128d y2v = cloud->gety2_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, x2v);
		_mm_storeh_pd(&x2, x2v);
		_mm_storel_pd(&y1, y2v);
		_mm_storeh_pd(&y2, y2v);

		const double dispX = x1 - x2;
		const double dispY = y1 - y2;
		const double disp = sqrt(dispX*dispX + dispY*dispY);
		__m128d dispXv = _mm_set_pd(dispX, dispX);
		__m128d dispYv = _mm_set_pd(dispY, dispY);
		__m128d dispV = _mm_set_pd(disp, disp);

		potential(currentParticle, disp, cloud->getq2r_pd(currentParticle));
		field(currentParticle, dispV, dispXv, dispYv, cloud->getphi2_pd(currentParticle));

		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			dispXv = x2v - cloud->getx2_pd(i);
			dispYv = y2v - cloud->gety2_pd(i);
			dispV = _mm_sqrt_pd(dispXv*dispXv + dispYv*dispYv);
			const __m128d dispXvr = x2v - cloud->getx2r_pd(i);
			const __m128d dispYvr = y2v - cloud->gety2r_pd(i);
			const __m128d dispVr = _mm_sqrt_pd(dispXvr*dispXvr + dispYvr*dispYvr);

			potential(currentParticle, i, dispV, cloud->getq2_pd(currentParticle), cloud->getq2_pd(i));
			potentialr(currentParticle, i, dispVr, cloud->getq2_pd(currentParticle), cloud->getq2r_pd(i));

			field(currentParticle, i, dispV, dispXv, dispYv, cloud->getphi2_pd(currentParticle), cloud->getphi2_pd(i));
			fieldr(currentParticle, i, dispVr, dispXvr, dispYvr, cloud->getphi2_pd(currentParticle), cloud->getphi2r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d x3v = cloud->getx3_pd(currentParticle);
		const __m128d y3v = cloud->gety3_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, x3v);
		_mm_storeh_pd(&x2, x3v);
		_mm_storel_pd(&y1, y3v);
		_mm_storeh_pd(&y2, y3v);

		const double dispX = x1 - x2;
		const double dispY = y1 - y2;
		const double disp = sqrt(dispX*dispX + dispY*dispY);
		__m128d dispXv = _mm_set_pd(dispX, dispX);
		__m128d dispYv = _mm_set_pd(dispY, dispY);
		__m128d dispV = _mm_set_pd(disp, disp);

		potential(currentParticle, disp, cloud->getq3r_pd(currentParticle));
		field(currentParticle, dispV, dispXv, dispYv, cloud->getphi3_pd(currentParticle));

		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			dispXv = x3v - cloud->getx3_pd(i);
			dispYv = y3v - cloud->gety3_pd(i);
			dispV = _mm_sqrt_pd(dispXv*dispXv + dispYv*dispYv);
			const __m128d dispXvr = x3v - cloud->getx3r_pd(i);
			const __m128d dispYvr = y3v - cloud->gety3r_pd(i);
			const __m128d dispVr = _mm_sqrt_pd(dispXvr*dispXvr + dispYvr*dispYvr);

			potential(currentParticle, i, dispV, cloud->getq3_pd(currentParticle), cloud->getq3_pd(i));
			potentialr(currentParticle, i, dispVr, cloud->getq3_pd(currentParticle), cloud->getq3r_pd(i));

			field(currentParticle, i, dispV, dispXv, dispYv, cloud->getphi3_pd(currentParticle), cloud->getphi3_pd(i));
			fieldr(currentParticle, i, dispVr, dispXvr, dispYvr, cloud->getphi3_pd(currentParticle), cloud->getphi3r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d x4v = cloud->getx4_pd(currentParticle);
		const __m128d y4v = cloud->gety4_pd(currentParticle);
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, x4v);
		_mm_storeh_pd(&x2, x4v);
		_mm_storel_pd(&y1, y4v);
		_mm_storeh_pd(&y2, y4v);

		const double dispX = x1 - x2;
		const double dispY = y1 - y2;
		const double disp = sqrt(dispX*dispX + dispY*dispY);
		__m128d dispXv = _mm_set_pd(dispX, dispX);
		__m128d dispYv = _mm_set_pd(dispY, dispY);
		__m128d dispV = _mm_set_pd(disp, disp);

		potential(currentParticle, disp, cloud->getq4r_pd(currentParticle));
		field(currentParticle, dispV, dispXv, dispYv, cloud->getphi4_pd(currentParticle));

		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			dispXv = x4v - cloud->getx4_pd(i);
			dispYv = y4v - cloud->gety4_pd(i);
			dispV = _mm_sqrt_pd(dispXv*dispXv + dispYv*dispYv);
			const __m128d dispXvr = x4v - cloud->getx4r_pd(i);
			const __m128d dispYvr = y4v - cloud->gety4r_pd(i);
			const __m128d dispVr = _mm_sqrt_pd(dispXvr*dispXvr + dispYvr*dispYvr);

			potential(currentParticle, i, dispV, cloud->getq4_pd(currentParticle), cloud->getq4_pd(i));
			potentialr(currentParticle, i, dispVr, cloud->getq4_pd(currentParticle), cloud->getq4r_pd(i));

			field(currentParticle, i, dispV, dispXv, dispYv, cloud->getphi4_pd(currentParticle), cloud->getphi4_pd(i));
			fieldr(currentParticle, i, dispVr, dispXvr, dispYvr, cloud->getphi4_pd(currentParticle), cloud->getphi4r_pd(i));
		}
	}
}

inline void FieldPotentialOperator::potential(const cloud_index currentParticle, const double disp, const __m128d charges)
{
	const double valExp = disp/cloud->dustDebye;

	if (valExp < 10.0) // restrict to 10*(ion debye length)
	{
		const double temp = 1.0/(4.0*M_PI*8.85E-12*disp*exp(valExp));

		cloud->phiCache[currentParticle] -= _mm_set_pd(temp, temp)*charges;
	}
}

inline void FieldPotentialOperator::potential(const cloud_index currentParticle, const cloud_index iParticle, const __m128d disp, const __m128d currentCharges, const __m128d iCharges)
{
	const __m128d valExp = disp/_mm_set_pd(cloud->dustDebye, cloud->dustDebye);
	
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
	
	__m128d expv = _mm_set_pd(boolH ? exp(expH) : 0.0, // _mm_set_pd is backwards
							  boolL ? exp(expL) : 0.0);

	// conclude potential calculation:
	const double k = 1.0/4.0*M_PI*8.85E-12;
	__m128d temp = _mm_set_pd(k, k)/(disp*expv);

	cloud->phiCache[currentParticle] -= temp*iCharges;
	cloud->phiCache[iParticle] -= temp*currentCharges;
}

inline void FieldPotentialOperator::potentialr(const cloud_index currentParticle, const cloud_index iParticle, const __m128d disp, const __m128d currentCharges, const __m128d iCharges)
{
	const __m128d valExp = disp/_mm_set_pd(cloud->dustDebye, cloud->dustDebye);
	
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
	
	__m128d expv = _mm_set_pd(boolH ? exp(expH) : 0.0, // _mm_set_pd is backwards
							  boolL ? exp(expL) : 0.0);

	// conclude potential calculation:
	const double k = 1.0/4.0*M_PI*8.85E-12;
	const __m128d temp = _mm_set_pd(k, k)/(disp*expv);
	const __m128d temp2 = temp*currentCharges;

	cloud->phiCache[currentParticle] -= temp*iCharges;
	cloud->phiCache[iParticle] -= _mm_shuffle_pd(temp2, temp2, _MM_SHUFFLE2(0, 1));
}

inline void FieldPotentialOperator::field(const cloud_index currentParticle, const __m128d dispV, const __m128d dispXv, const __m128d dispYv, const __m128d phi)
{
		const double debye = cloud->dustDebye;
		const __m128d temp = (_mm_set_pd(1.0,1.0) + dispV/_mm_set_pd(debye, debye))/(dispV*dispV)*phi;
		cloud->ExCache[currentParticle] += temp*dispXv;
		cloud->EyCache[currentParticle] += temp*dispYv;
}

inline void FieldPotentialOperator::field(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispV, const __m128d dispXv, const __m128d dispYv, const __m128d currentPhi, const __m128d iPhi)
{
	const double debye = cloud->dustDebye;
	const __m128d temp = (_mm_set_pd(1.0,1.0) + dispV/_mm_set_pd(debye, debye))/(dispV*dispV);

	cloud->ExCache[currentParticle] += temp*currentPhi*dispXv;
	cloud->EyCache[currentParticle] += temp*currentPhi*dispYv;

	cloud->ExCache[iParticle] += temp*iPhi*dispXv;
	cloud->EyCache[iParticle] += temp*iPhi*dispYv;
}

inline void FieldPotentialOperator::fieldr(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispVr, const __m128d dispXvr, const __m128d dispYvr, const __m128d currentPhi, const __m128d iPhir)
{
	const double debye = cloud->dustDebye;
	const __m128d temp = (_mm_set_pd(1.0,1.0) + dispVr/_mm_set_pd(debye, debye))/(dispVr*dispVr);
	const __m128d temp2 = temp*iPhir;
	const __m128d temp3 = temp2*dispXvr;
	const __m128d temp4 = temp2*dispYvr;

	cloud->ExCache[currentParticle] += temp*currentPhi*dispXvr;
	cloud->EyCache[currentParticle] += temp*currentPhi*dispYvr;

	cloud->ExCache[iParticle] += _mm_shuffle_pd(temp3, temp3, _MM_SHUFFLE2(0, 1));
	cloud->EyCache[iParticle] += _mm_shuffle_pd(temp4, temp4, _MM_SHUFFLE2(0, 1));
}
