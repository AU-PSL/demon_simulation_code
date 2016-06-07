/**
* @file  CacheOperator.cpp
* @class CacheOperator CacheOperator.h
*
* @brief Caches positions of velocities for the Runge-Kutta substeps
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "CacheOperator.h"

/**
* @brief The first Runge-Kutta sub-timeStep positions and velocities remain unaltered.
*
* @param[in] currentTime Current simulaton time
**/
void CacheOperator::operation1(const double currentTime) {
    (void)currentTime;
}

/**
* @brief Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
*        take the form of (a + b1/2)
*
* @param[in] currentTime Current simulaton time
**/
void CacheOperator::operation2(const double currentTime) {
    (void)currentTime;
	const doubleV halfv = set1_pd(0.5);
    BEGIN_PARALLEL_FOR(i, e, cloud->n/DOUBLE_STRIDE, 1, static)
		const cloud_index offset = DOUBLE_STRIDE*i;
		cloud->xCache[i] = fmadd_pd(halfv, load_pd(cloud->l1 + offset), load_pd(cloud->x + offset));
		cloud->yCache[i] = fmadd_pd(halfv, load_pd(cloud->n1 + offset), load_pd(cloud->y + offset));
		cloud->VxCache[i] = fmadd_pd(halfv, load_pd(cloud->k1 + offset), load_pd(cloud->Vx + offset));
		cloud->VyCache[i] = fmadd_pd(halfv, load_pd(cloud->m1 + offset), load_pd(cloud->Vy + offset));
    END_PARALLEL_FOR
}

/**
* @brief Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
*        take the form of (a + b2/2)
*
* @param[in] currentTime Current simulaton time
**/
void CacheOperator::operation3(const double currentTime) {
    (void)currentTime;
	const doubleV halfv = set1_pd(0.5);
	BEGIN_PARALLEL_FOR(i, e, cloud->n/DOUBLE_STRIDE, 1, static)
		const cloud_index offset = DOUBLE_STRIDE*i;
		cloud->xCache[i] = fmadd_pd(halfv, load_pd(cloud->l2 + offset), load_pd(cloud->x + offset));
		cloud->yCache[i] = fmadd_pd(halfv, load_pd(cloud->n2 + offset), load_pd(cloud->y + offset));
		cloud->VxCache[i] = fmadd_pd(halfv, load_pd(cloud->k2 + offset), load_pd(cloud->Vx + offset));
		cloud->VyCache[i] = fmadd_pd(halfv, load_pd(cloud->m2 + offset), load_pd(cloud->Vy + offset));
	END_PARALLEL_FOR
}

/**
* @brief Cache positions and velocities for second Runge-Kutta sub-timeStep. Values 
*        take the form of (a + b3)
*
* @param[in] currentTime Current simulaton time
**/
void CacheOperator::operation4(const double currentTime) {
    (void)currentTime;
	BEGIN_PARALLEL_FOR(i, e, cloud->n/DOUBLE_STRIDE, 1, static)
		const cloud_index offset = DOUBLE_STRIDE*i;
		cloud->xCache[i] = add_pd(load_pd(cloud->x + offset), load_pd(cloud->l3 + offset));
		cloud->yCache[i] = add_pd(load_pd(cloud->y + offset), load_pd(cloud->n3 + offset));
		cloud->VxCache[i] = add_pd(load_pd(cloud->Vx + offset), load_pd(cloud->k3 + offset));
		cloud->VyCache[i] = add_pd(load_pd(cloud->Vy + offset), load_pd(cloud->m3 + offset));
	END_PARALLEL_FOR
}
