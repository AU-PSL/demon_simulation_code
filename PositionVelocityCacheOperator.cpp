/*===- PositionVelocityCacheOperator.cpp - libSimulation -======================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "PositionVelocityCacheOperator.h"
#include "VectorCompatibility.h"

void PositionVelocityCacheOperator::operation1(const double currentTime) 
{
	// For the first RK4 timeStep the position and velocity remain unaltered.
}

void PositionVelocityCacheOperator::operation2(const double currentTime) 
{
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l1 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n1 + offset)/twov;
		cloud->zCache[i] = _mm_load_pd(cloud->z + offset) + _mm_load_pd(cloud->p1 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k1 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m1 + offset)/twov;
		cloud->VzCache[i] = _mm_load_pd(cloud->Vz + offset) + _mm_load_pd(cloud->o1 + offset)/twov;
	}
}

void PositionVelocityCacheOperator::operation3(const double currentTime) 
{
	const __m128d twov = _mm_set1_pd(2.0);
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l2 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n2 + offset)/twov;
		cloud->zCache[i] = _mm_load_pd(cloud->z + offset) + _mm_load_pd(cloud->p2 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k2 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m2 + offset)/twov;
		cloud->VzCache[i] = _mm_load_pd(cloud->Vz + offset) + _mm_load_pd(cloud->o2 + offset)/twov;
	}
}

void PositionVelocityCacheOperator::operation4(const double currentTime) 
{
	for (cloud_index i = 0, e = cloud->n/2; i < e; i++) 
	{
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l3 + offset);
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n3 + offset);
		cloud->zCache[i] = _mm_load_pd(cloud->z + offset) + _mm_load_pd(cloud->p3 + offset);
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k3 + offset);
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m3 + offset);
		cloud->VzCache[i] = _mm_load_pd(cloud->Vz + offset) + _mm_load_pd(cloud->o3 + offset);
	}
}