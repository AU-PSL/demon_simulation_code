/*===- Runge_Kutta2.cpp - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "Runge_Kutta2.h"

Runge_Kutta2::Runge_Kutta2(Cloud * const C, const ForceArray &FA, 
                           const double timeStep, const double startTime)
: Integrator(C, FA, timeStep, startTime) {}

// 4th order Runge-Kutta algorithm:
void Runge_Kutta2::moveParticles(const double endTime) {
	while (currentTime < endTime) {
		// Second argument must be 2 more than the first.
		const double dt = modifyTimeStep(1.0e-4f, init_dt); // implement dynamic timstep (if necessary):
		const __m128d vdt = _mm_set1_pd(dt); // store timestep as vector const
		
		const cloud_index numParticles = cloud->n;
        
		operate1(currentTime);
		force1(currentTime); // compute net force1
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2, static) // calculate k1 and l1 for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass into vector
            
			// assign force pointers for stylistic purposes:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
            
			// calculate ith and (i+1)th tidbits: 
			_mm_store_pd(cloud->k1 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l1 + i, vdt*cloud->getVx1_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m1 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n1 + i, vdt*cloud->getVy1_pd(i)); // positionY tidbit
            
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
		END_PARALLEL_FOR
        
		operate2(currentTime + dt/2.0);
		force2(currentTime + dt/2.0); // compute net force2
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2, static) // calculate k2 and l for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass
            
			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
            
			// calculate ith and (i+1)th tidbits: 
			_mm_store_pd(cloud->k2 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l2 + i, vdt*cloud->getVx2_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m2 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n2 + i, vdt*cloud->getVy2_pd(i)); // positionY tidbit
            
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
		END_PARALLEL_FOR
        
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2, static) // calculate next position and next velocity for entire cloud
			// load ith and (i+1)th k's into vectors:
			const __m128d vk2 = _mm_load_pd(cloud->k2 + i);
			const __m128d vl2 = _mm_load_pd(cloud->l2 + i);
			const __m128d vm2 = _mm_load_pd(cloud->m2 + i);
			const __m128d vn2 = _mm_load_pd(cloud->n2 + i);
            
			// assign position and velocity pointers (stylistic):
			double * const px = cloud->x + i;
			double * const py = cloud->y + i;
			double * const pVx = cloud->Vx + i;
			double * const pVy = cloud->Vy + i;
            
			// calculate next positions and velocities:
			_mm_store_pd(pVx, _mm_load_pd(pVx) + vk2);
			_mm_store_pd(px, _mm_load_pd(px) + vl2);
			_mm_store_pd(pVy, _mm_load_pd(pVy) + vm2);
			_mm_store_pd(py, _mm_load_pd(py) + vn2);
		END_PARALLEL_FOR
        
		currentTime += dt;
	}
}

void Runge_Kutta2::operate1(const double time) const {
 	for (Operator * const O : operations)
		O->operation1(time);
}

void Runge_Kutta2::operate2(const double time) const {
 	for (Operator * const O : operations)
		O->operation2(time);
}

void Runge_Kutta2::force1(const double time) const {
	for (Force * const F : forces)
		F->force1(time);
}

void Runge_Kutta2::force2(const double time) const {
	for (Force * const F : forces)
		F->force2(time);
}
