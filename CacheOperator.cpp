/*===- CacheOperator.cpp - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "CacheOperator.h"
#include "VectorCompatibility.h"

void CacheOperator::operation1(const double currentTime) 
{
    (void)currentTime;
	// For the first RK4 timeStep the position and velocity remain unaltered.
}

void CacheOperator::operation2(const double currentTime) 
{
    (void)currentTime;
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l1 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n1 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k1 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m1 + offset)/twov;
		
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q1 + offset)/twov;
	}
}

void CacheOperator::operation3(const double currentTime) 
{
    (void)currentTime;
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l2 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n2 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k2 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m2 + offset)/twov;
		
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q2 + offset)/twov;
	}
}

void CacheOperator::operation4(const double currentTime) 
{
    (void)currentTime;
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l3 + offset);
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n3 + offset);
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k3 + offset);
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m3 + offset);
		
		cloud->qCache[i] = _mm_load_pd(cloud->charge + offset) + _mm_load_pd(cloud->q3 + offset);
	}
}
