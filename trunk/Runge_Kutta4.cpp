/*===- Runge_Kutta4.cpp - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "Runge_Kutta4.h"

Runge_Kutta4::Runge_Kutta4(Cloud * const C, const ForceArray &FA, 
                           const double timeStep, const double startTime)
: Runge_Kutta2(C, FA, timeStep, startTime) {}

// 4th order Runge-Kutta algorithm:
void Runge_Kutta4::moveParticles(const double endTime) {
	// create vector constants:
	const doubleV v2 = _mm_set1_pd(2.0);
	const doubleV v6 = _mm_set1_pd(6.0);
    
	while (currentTime < endTime) {
		const double dt = modifyTimeStep(1.0e-4f, init_dt); // implement dynamic timstep (if necessary):
		const doubleV vdt = _mm_set1_pd(dt); // store timestep as vector const
		
		const cloud_index numParticles = cloud->n;
        
		operate1(currentTime);
		force1(currentTime); // compute net force1
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k1 and l1 for entire cloud
			const doubleV vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass into vector

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
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k2 and l2 for entire cloud
			const doubleV vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

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

		operate3(currentTime + dt/2.0);
		force3(currentTime + dt/2.0); // compute net force3
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k3 and l3 for entire cloud
			const doubleV vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;

			// calculate ith and (i+1)th tibits: 
			_mm_store_pd(cloud->k3 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l3 + i, vdt*cloud->getVx3_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m3 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n3 + i, vdt*cloud->getVy3_pd(i)); // positionY tidbit
			
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
		END_PARALLEL_FOR
        
		operate4(currentTime + dt);
		force4(currentTime + dt); // compute net force4
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k4 and l4 for entire cloud
			const doubleV vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
            
			_mm_store_pd(cloud->k4 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l4 + i, vdt*cloud->getVx4_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m4 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n4 + i, vdt*cloud->getVy4_pd(i)); // positionY tidbit
			
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
		END_PARALLEL_FOR

        // Calculate next position and next velocity for entire cloud.
        BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static)
			// load ith and (i+1)th k's into vectors:
			const doubleV vk1 = _mm_load_pd(cloud->k1 + i);
			const doubleV vk2 = _mm_load_pd(cloud->k2 + i);
			const doubleV vk3 = _mm_load_pd(cloud->k3 + i);
			const doubleV vk4 = _mm_load_pd(cloud->k4 + i);

			// load ith and (i+1)th l's into vectors: 
			const doubleV vl1 = _mm_load_pd(cloud->l1 + i);
			const doubleV vl2 = _mm_load_pd(cloud->l2 + i);
			const doubleV vl3 = _mm_load_pd(cloud->l3 + i);
			const doubleV vl4 = _mm_load_pd(cloud->l4 + i);

			// load ith and (i+1)th m's into vectors: 
			const doubleV vm1 = _mm_load_pd(cloud->m1 + i);
			const doubleV vm2 = _mm_load_pd(cloud->m2 + i);
			const doubleV vm3 = _mm_load_pd(cloud->m3 + i);
			const doubleV vm4 = _mm_load_pd(cloud->m4 + i);

			// load ith and (i+1)th n's into vectors:
			const doubleV vn1 = _mm_load_pd(cloud->n1 + i);
			const doubleV vn2 = _mm_load_pd(cloud->n2 + i);
			const doubleV vn3 = _mm_load_pd(cloud->n3 + i);
			const doubleV vn4 = _mm_load_pd(cloud->n4 + i);

			// calculate next positions and velocities:
			plusEqual_pd(cloud->Vx + i, (vk1 + v2*(vk2 + vk3) + vk4)/v6);
			plusEqual_pd(cloud->x + i, (vl1 + v2*(vl2 + vl3) + vl4)/v6);
			plusEqual_pd(cloud->Vy + i, (vm1 + v2*(vm2 + vm3) + vm4)/v6);
			plusEqual_pd(cloud->y + i, (vn1 + v2*(vn2 + vn3) + vn4)/v6);
		END_PARALLEL_FOR

		currentTime += dt;
	}
}

inline void Runge_Kutta4::operate3(const double time) const {
 	for (Operator * const O : operations)
		O->operation3(time);
}

inline void Runge_Kutta4::operate4(const double time) const {
 	for (Operator * const O : operations)
		O->operation4(time);
}

inline void Runge_Kutta4::force3(const double time) const {
	for (Force * const F : forces)
		F->force3(time);
}

inline void Runge_Kutta4::force4(const double time) const {
	for (Force * const F : forces)
		F->force4(time);
}
