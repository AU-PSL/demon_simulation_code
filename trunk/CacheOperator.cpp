/*===- CacheOperator.cpp - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "CacheOperator.h"

// The first Runge-Kutta sub-timeStep positions and velocitys remain unaltered.
void CacheOperator::operation1(const double currentTime) {
    (void)currentTime;
}

// Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
// take the form of (a + b1/2)
void CacheOperator::operation2(const double currentTime) {
    (void)currentTime;
	const doubleV halfv = _mm_set1_pd(0.5);
    BEGIN_PARALLEL_FOR(i, e, cloud->n/2, 1, static)
		const cloud_index offset = 2*i;
		cloud->xCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->l1 + offset), _mm_load_pd(cloud->x + offset));
		cloud->yCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->n1 + offset), _mm_load_pd(cloud->y + offset));
		cloud->VxCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->k1 + offset), _mm_load_pd(cloud->Vx + offset));
		cloud->VyCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->m1 + offset), _mm_load_pd(cloud->Vy + offset));
    END_PARALLEL_FOR
}

// Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
// take the form of (a + b2/2)
void CacheOperator::operation3(const double currentTime) {
    (void)currentTime;
	const doubleV halfv = _mm_set1_pd(0.5);
	BEGIN_PARALLEL_FOR(i, e, cloud->n/2, 1, static)
		const cloud_index offset = 2*i;
		cloud->xCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->l2 + offset), _mm_load_pd(cloud->x + offset));
		cloud->yCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->n2 + offset), _mm_load_pd(cloud->y + offset));
		cloud->VxCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->k2 + offset), _mm_load_pd(cloud->Vx + offset));
		cloud->VyCache[i] = fmadd_pd(halfv, _mm_load_pd(cloud->m2 + offset), _mm_load_pd(cloud->Vy + offset));
	END_PARALLEL_FOR
}

// Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
// take the form of (a + b3)
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
