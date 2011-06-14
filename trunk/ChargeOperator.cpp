/*===- ChargeOperator.cpp - libSimulation -=====================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "ChargeOperator.h"
#include "VectorCompatibility.h"

void ChargeOperator::operation1(const double currentTime) 
{
	// For the first RK4 timestep the charges remain unaltered.
}

void ChargeOperator::operation2(const double currentTime) 
{
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
//		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q1 + offset)/twov;
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset);
	}
}

void ChargeOperator::operation3(const double currentTime) 
{
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
//		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q2 + offset)/twov;
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset);
	}
}

void ChargeOperator::operation4(const double currentTime) 
{
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
//		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q3 + offset);
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset);
	}
}
