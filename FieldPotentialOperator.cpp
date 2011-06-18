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

const __m128d FieldPotentialOperator::coulomb = _mm_set1_pd(1.0/(4.0*M_PI*8.8542E-12));

void FieldPotentialOperator::operation1(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx1_pd(currentParticle);
		const __m128d vy1 = cloud->gety1_pd(currentParticle);
		const __m128d vq1 = cloud->getq1_pd(currentParticle);
		
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		
		fieldsAndPotentials(currentParticle/2, cloud->getq1r_pd(currentParticle), x1 - x2, y1 - y2);
		
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			fieldsAndPotentials(currentParticle/2, i/2, vq1, cloud->getq1_pd(i), vx1 - cloud->getx1_pd(i), vy1 - cloud->gety1_pd(i));
			fieldsAndPotentialsr(currentParticle/2, i/2, vq1, cloud->getq1r_pd(i), vx1 - cloud->getx1r_pd(i), vy1 - cloud->gety1r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation2(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx2_pd(currentParticle);
		const __m128d vy1 = cloud->gety2_pd(currentParticle);
		const __m128d vq1 = cloud->getq2_pd(currentParticle);
		
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		
		fieldsAndPotentials(currentParticle/2, cloud->getq2r_pd(currentParticle), x1 - x2, y1 - y2);
		
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			fieldsAndPotentials(currentParticle/2, i/2, vq1, cloud->getq1_pd(i), vx1 - cloud->getx2_pd(i), vy1 - cloud->gety2_pd(i));
			fieldsAndPotentialsr(currentParticle/2, i/2, vq1, cloud->getq1r_pd(i), vx1 - cloud->getx2r_pd(i), vy1 - cloud->gety2r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation3(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx3_pd(currentParticle);
		const __m128d vy1 = cloud->gety3_pd(currentParticle);
		const __m128d vq1 = cloud->getq3_pd(currentParticle);
		
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		
		fieldsAndPotentials(currentParticle/2, cloud->getq3r_pd(currentParticle), x1 - x2, y1 - y2);
		
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			fieldsAndPotentials(currentParticle/2, i/2, vq1, cloud->getq3_pd(i), vx1 - cloud->getx3_pd(i), vy1 - cloud->gety3_pd(i));
			fieldsAndPotentialsr(currentParticle/2, i/2, vq1, cloud->getq3r_pd(i), vx1 - cloud->getx3r_pd(i), vy1 - cloud->gety3r_pd(i));
		}
	}
}

void FieldPotentialOperator::operation4(const double currentTime)
{
	for (cloud_index currentParticle = 0, numParticles = cloud->n, e = cloud->n - 1; currentParticle < e; currentParticle += 2) 
	{
		const __m128d vx1 = cloud->getx4_pd(currentParticle);
		const __m128d vy1 = cloud->gety4_pd(currentParticle);
		const __m128d vq1 = cloud->getq4_pd(currentParticle);
		
		double x1, x2, y1, y2;
		_mm_storel_pd(&x1, vx1);
		_mm_storeh_pd(&x2, vx1);
		_mm_storel_pd(&y1, vy1);
		_mm_storeh_pd(&y2, vy1);
		
		fieldsAndPotentials(currentParticle/2, cloud->getq4r_pd(currentParticle), x1 - x2, y1 - y2);
		
		for (cloud_index i = currentParticle + 2; i < numParticles; i += 2)
		{
			fieldsAndPotentials(currentParticle/2, i/2, vq1, cloud->getq4_pd(i), vx1 - cloud->getx4_pd(i), vy1 - cloud->gety4_pd(i));
			fieldsAndPotentialsr(currentParticle/2, i/2, vq1, cloud->getq4r_pd(i), vx1 - cloud->getx4r_pd(i), vy1 - cloud->gety4r_pd(i));
		}
	}
}

inline void FieldPotentialOperator::fieldsAndPotentials(const cloud_index currentParticle, const __m128d currentCharge,
                                                        const double displacementX, const double displacementY)
{
	const double displacement = sqrt(displacementX*displacementX + displacementY*displacementY);
	const double valExp = displacement/cloud->dustDebye;
	
	if (valExp < 10.0) // restrict to 10*(ion debye length)
	{
		const __m128d coeffient = coulomb/_mm_set1_pd(exp(valExp));
		
		cloud->phi[currentParticle] += currentCharge/_mm_set1_pd(displacement)*coeffient;
		
		const __m128d eField = currentCharge*_mm_set1_pd((1 + displacement/cloud->dustDebye)/(displacement*displacement*displacement))*coeffient;
		cloud->Ex[currentParticle] += eField*_mm_set_pd(-displacementX, displacementX);
		cloud->Ey[currentParticle] += eField*_mm_set_pd(-displacementY, displacementY);
	}
}

inline void FieldPotentialOperator::fieldsAndPotentials(const cloud_index currentParticle, const cloud_index iParticle,
                                                        const __m128d currentCharge, const __m128d iCharge,
                                                        const __m128d displacementX, const __m128d displacementY)
{
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement/_mm_set1_pd(cloud->dustDebye);
	
	double valExpL, valExpH;
	_mm_storel_pd(&valExpL, valExp);
	_mm_storeh_pd(&valExpH, valExp);
	
	const bool boolL = valExpL < 10.0;
	const bool boolH = valExpH < 10.0;
	if (!boolL && !boolH)
		return;
	
	const __m128d expv = _mm_set_pd(boolH ? exp(-valExpH) : 0.0,
									boolL ? exp(-valExpL) : 0.0);
	const __m128d coeffient = coulomb*expv;
	
	cloud->phi[currentParticle] += iCharge/displacement*coeffient;
	cloud->phi[iParticle] += currentCharge/displacement*coeffient;
	
	const __m128d coeffient2 = coeffient*(_mm_set1_pd(1) + displacement/_mm_set1_pd(cloud->dustDebye))/(displacement*displacement*displacement);
	const __m128d xunit = displacementX*coeffient2;
	const __m128d yunit = displacementY*coeffient2;
	cloud->Ex[currentParticle] += iCharge*xunit;
	cloud->Ey[currentParticle] += iCharge*yunit;
	cloud->Ex[iParticle] -= currentCharge*xunit;
	cloud->Ey[iParticle] -= currentCharge*yunit;
}

inline void FieldPotentialOperator::fieldsAndPotentialsr(const cloud_index currentParticle, const cloud_index iParticle,
                                                         const __m128d currentCharge, const __m128d iCharge,
                                                         const __m128d displacementX, const __m128d displacementY)
{
	const __m128d displacement = _mm_sqrt_pd(displacementX*displacementX + displacementY*displacementY);
	const __m128d valExp = displacement/_mm_set1_pd(cloud->dustDebye);
	
	double valExpL, valExpH;
	_mm_storel_pd(&valExpL, valExp);
	_mm_storeh_pd(&valExpH, valExp);
	
	const bool boolL = valExpL < 10.0;
	const bool boolH = valExpH < 10.0;
	if (!boolL && !boolH)
		return;
	
	const __m128d expv = _mm_set_pd(boolH ? exp(-valExpH) : 0.0,
									boolL ? exp(-valExpL) : 0.0);
	const __m128d coeffient = coulomb*expv;
	
	cloud->phi[currentParticle] += iCharge/displacement*coeffient;
	const __m128d iPhi = currentCharge/displacement*coeffient;
	cloud->phi[iParticle] += _mm_shuffle_pd(iPhi, iPhi, _MM_SHUFFLE2(0, 1));
	
	const __m128d coeffient2 = coeffient*(_mm_set1_pd(1) + displacement/_mm_set1_pd(cloud->dustDebye))/(displacement*displacement*displacement);
	const __m128d xunit = displacementX*coeffient2;
	const __m128d yunit = displacementY*coeffient2;
	cloud->Ex[currentParticle] += iCharge*xunit;
	cloud->Ey[currentParticle] += iCharge*yunit;
	
	const __m128d iEx = currentCharge*xunit;
	const __m128d iEy = currentCharge*yunit;
	cloud->Ex[iParticle] -= _mm_shuffle_pd(iEx, iEx, _MM_SHUFFLE2(0, 1));
	cloud->Ey[iParticle] -= _mm_shuffle_pd(iEy, iEy, _MM_SHUFFLE2(0, 1));
}
