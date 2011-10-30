/*===- CacheOperator.cpp - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "CacheOperator.h"

void CacheOperator::operation1(const double currentTime) {
    (void)currentTime;
	// For the first RK4 timeStep the position and velocity remain unaltered.
}

void CacheOperator::operation2(const double currentTime) {
    (void)currentTime;
	const __m128d twov = _mm_set1_pd(2.0);
    BEGIN_PARALLEL_FOR(i, e, cloud->n/2, 1, static)
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l1 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n1 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k1 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m1 + offset)/twov;
    END_PARALLEL_FOR
}

void CacheOperator::operation3(const double currentTime) {
    (void)currentTime;
	const __m128d twov = _mm_set1_pd(2.0);
	BEGIN_PARALLEL_FOR(i, e, cloud->n/2, 1, static)
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l2 + offset)/twov;
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n2 + offset)/twov;
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k2 + offset)/twov;
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m2 + offset)/twov;
	END_PARALLEL_FOR
}

void CacheOperator::operation4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(i, e, cloud->n/2, 1, static)
		const cloud_index offset = 2*i;
		cloud->xCache[i] = _mm_load_pd(cloud->x + offset) + _mm_load_pd(cloud->l3 + offset);
		cloud->yCache[i] = _mm_load_pd(cloud->y + offset) + _mm_load_pd(cloud->n3 + offset);
		cloud->VxCache[i] = _mm_load_pd(cloud->Vx + offset) + _mm_load_pd(cloud->k3 + offset);
		cloud->VyCache[i] = _mm_load_pd(cloud->Vy + offset) + _mm_load_pd(cloud->m3 + offset);
	END_PARALLEL_FOR
}
