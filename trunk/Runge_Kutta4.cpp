/*===- Runge_Kutta4.cpp - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "Runge_Kutta4.h"

Runge_Kutta4::Runge_Kutta4(Cloud * const myCloud, Force **forces, const force_index forcesSize, 
                           const double timeStep, const double startTime)
: Runge_Kutta2(myCloud, forces, forcesSize, timeStep, startTime) {}

// 4th order Runge-Kutta algorithm:
void Runge_Kutta4::moveParticles(const double endTime) {
	// create vector constants:
	const __m128d v2 = _mm_set1_pd(2.0);
	const __m128d v6 = _mm_set1_pd(6.0);

#ifdef CHARGE
	const __m128d qConst3 = _mm_set1_pd(4.0*M_PI*Cloud::particleRadius*Cloud::epsilon0);
#endif
    
	while (currentTime < endTime) {
		// Second argument must be 2 more than the first.
		const double dt = modifyTimeStep(1.0e-4, init_dt); // implement dynamic timstep (if necessary):
		const __m128d vdt = _mm_set1_pd(dt); // store timestep as vector const
		
		const cloud_index numParticles = cloud->n;
        
		operate1(currentTime);
		force1(currentTime); // compute net force1
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2) // calculate k1 and l1 for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass into vector

			// assign force pointers for stylistic purposes:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
			double * const pPhi = cloud->phi + i;
           
			// calculate ith and (i+1)th tidbits: 
			_mm_store_pd(cloud->k1 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l1 + i, vdt*cloud->getVx1_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m1 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n1 + i, vdt*cloud->getVy1_pd(i)); // positionY tidbit
#ifdef CHARGE
			const __m128d pQ = cloud->getq1_pd(i);
            __m128d qConst1, qConst2;
            Cloud::setChargeConsts(pQ, qConst1, qConst2);
			_mm_store_pd(cloud->q1 + i, -vdt*(qConst1*pQ + qConst2*qConst3*_mm_load_pd(pPhi)));
#else
			_mm_store_pd(cloud->q1 + i, _mm_setzero_pd()); // charge tidbit
#endif

			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
			_mm_store_pd(pPhi, _mm_setzero_pd());
		END_PARALLEL_FOR
        
		operate2(currentTime + dt/2.0);
		force2(currentTime + dt/2.0); // compute net force2
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2) // calculate k2 and l2 for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
			double * const pPhi = cloud->phi + i;

			// calculate ith and (i+1)th tidbits: 
			_mm_store_pd(cloud->k2 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l2 + i, vdt*cloud->getVx2_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m2 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n2 + i, vdt*cloud->getVy2_pd(i)); // positionY tidbit
#ifdef CHARGE
			const __m128d pQ = cloud->getq2_pd(i);
			__m128d qConst1, qConst2;
			Cloud::setChargeConsts(pQ, qConst1, qConst2);
			_mm_store_pd(cloud->q2 + i, -vdt*(qConst1*pQ + qConst2*qConst3*_mm_load_pd(pPhi)));
#else			
			_mm_store_pd(cloud->q2 + i, _mm_setzero_pd()); // charge tidbit
#endif

			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
			_mm_store_pd(pPhi, _mm_setzero_pd());
		END_PARALLEL_FOR

		operate3(currentTime + dt/2.0);
		force3(currentTime + dt/2.0); // compute net force3
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2) // calculate k3 and l3 for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
			double * const pPhi = cloud->phi + i;

			// calculate ith and (i+1)th tibits: 
			_mm_store_pd(cloud->k3 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l3 + i, vdt*cloud->getVx3_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m3 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n3 + i, vdt*cloud->getVy3_pd(i)); // positionY tidbit
#ifdef CHARGE
			const __m128d pQ = cloud->getq3_pd(i);
			__m128d qConst1, qConst2;
			Cloud::setChargeConsts(pQ, qConst1, qConst2);
			_mm_store_pd(cloud->q3 + i, -vdt*(qConst1*pQ + qConst2*qConst3*_mm_load_pd(pPhi)));
#else			
			_mm_store_pd(cloud->q3 + i, _mm_setzero_pd()); // charge tidbit
#endif
			
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
			_mm_store_pd(pPhi, _mm_setzero_pd());
		END_PARALLEL_FOR
        
		operate4(currentTime + dt);
		force4(currentTime + dt); // compute net force4
		BEGIN_PARALLEL_FOR(i, e, numParticles, 2) // calculate k4 and l4 for entire cloud
			const __m128d vmass = _mm_load_pd(cloud->mass + i); // load ith and (i+1)th mass

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
			double * const pPhi = cloud->phi + i;
            
			_mm_store_pd(cloud->k4 + i, vdt*_mm_load_pd(pFx)/vmass); // velocityX tidbit
			_mm_store_pd(cloud->l4 + i, vdt*cloud->getVx4_pd(i)); // positionX tidbit
			_mm_store_pd(cloud->m4 + i, vdt*_mm_load_pd(pFy)/vmass); // velocityY tidbit
			_mm_store_pd(cloud->n4 + i, vdt*cloud->getVy4_pd(i)); // positionY tidbit
#ifdef CHARGE
			const __m128d pQ = cloud->getq4_pd(i);
			__m128d qConst1, qConst2;
			Cloud::setChargeConsts(pQ, qConst1, qConst2);
			_mm_store_pd(cloud->q4 + i, -vdt*(qConst1*pQ + qConst2*qConst3*_mm_load_pd(pPhi)));
#else			
			_mm_store_pd(cloud->q4 + i, _mm_setzero_pd()); // charge tidbit
#endif
			
			// reset forces to zero:
			_mm_store_pd(pFx, _mm_setzero_pd());
			_mm_store_pd(pFy, _mm_setzero_pd());
			_mm_store_pd(pPhi, _mm_setzero_pd());
		END_PARALLEL_FOR

        BEGIN_PARALLEL_FOR(i, e, numParticles, 2) // calculate next position and next velocity for entire cloud
			// load ith and (i+1)th k's into vectors:
			const __m128d vk1 = _mm_load_pd(cloud->k1 + i);
			const __m128d vk2 = _mm_load_pd(cloud->k2 + i);
			const __m128d vk3 = _mm_load_pd(cloud->k3 + i);
			const __m128d vk4 = _mm_load_pd(cloud->k4 + i);

			// load ith and (i+1)th l's into vectors: 
			const __m128d vl1 = _mm_load_pd(cloud->l1 + i);
			const __m128d vl2 = _mm_load_pd(cloud->l2 + i);
			const __m128d vl3 = _mm_load_pd(cloud->l3 + i);
			const __m128d vl4 = _mm_load_pd(cloud->l4 + i);

			// load ith and (i+1)th m's into vectors: 
			const __m128d vm1 = _mm_load_pd(cloud->m1 + i);
			const __m128d vm2 = _mm_load_pd(cloud->m2 + i);
			const __m128d vm3 = _mm_load_pd(cloud->m3 + i);
			const __m128d vm4 = _mm_load_pd(cloud->m4 + i);

			// load ith and (i+1)th n's into vectors:
			const __m128d vn1 = _mm_load_pd(cloud->n1 + i);
			const __m128d vn2 = _mm_load_pd(cloud->n2 + i);
			const __m128d vn3 = _mm_load_pd(cloud->n3 + i);
			const __m128d vn4 = _mm_load_pd(cloud->n4 + i);
			
			// load ith and (i+1)th q's into vectors:
			const __m128d vq1 = _mm_load_pd(cloud->q1 + i);
			const __m128d vq2 = _mm_load_pd(cloud->q2 + i);
			const __m128d vq3 = _mm_load_pd(cloud->q3 + i);
			const __m128d vq4 = _mm_load_pd(cloud->q4 + i);

			// assign position and velocity pointers (stylistic):
			double * const px = cloud->x + i;
			double * const py = cloud->y + i;
			double * const pVx = cloud->Vx + i;
			double * const pVy = cloud->Vy + i;
			
			double * const pC = cloud->charge + i;

			// calculate next positions and velocities:
			_mm_store_pd(pVx, _mm_load_pd(pVx) + (vk1 + v2*(vk2 + vk3) + vk4)/v6);
			_mm_store_pd(px, _mm_load_pd(px) + (vl1 + v2*(vl2 + vl3) + vl4)/v6);
			_mm_store_pd(pVy, _mm_load_pd(pVy) + (vm1 + v2*(vm2 + vm3) + vm4)/v6);
			_mm_store_pd(py, _mm_load_pd(py) + (vn1 + v2*(vn2 + vn3) + vn4)/v6);
			
			_mm_store_pd(pC, _mm_load_pd(pC) + (vq1 + v2*(vq2 + vq3) + vq4)/v6);
		END_PARALLEL_FOR

		currentTime += dt;
	}
}

inline void Runge_Kutta4::operate3(const double time) const {
 	for (operator_index i = 0; i < numOperators; i++)
		operations[i]->operation3(time);
}

inline void Runge_Kutta4::operate4(const double time) const {
 	for (operator_index i = 0; i < numOperators; i++)
		operations[i]->operation4(time);
}

inline void Runge_Kutta4::force3(const double time) const {
 	for (force_index i = 0; i < numForces; i++)
		theForce[i]->force3(time);
}

inline void Runge_Kutta4::force4(const double time) const {
 	for (force_index i = 0; i < numForces; i++)
		theForce[i]->force4(time);
}
