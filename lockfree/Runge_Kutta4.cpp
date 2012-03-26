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
	while (currentTime < endTime) {
		const double dt = modifyTimeStep(1.0e-4f, init_dt); // implement dynamic timstep (if necessary):
		const doubleV vdt = set1_pd(dt); // store timestep as vector const
		
		const cloud_index numParticles = cloud->n;
        
		operate1(currentTime);
		force1(currentTime); // compute net force1
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k1 and l1 for entire cloud
			const doubleV vmass = load_pd(cloud->mass + i);

			// assign force pointers for stylistic purposes:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
           
            store_pd(cloud->k1 + i, div_pd(mul_pd(vdt, load_pd(pFx)), vmass)); // velocityX tidbit
            store_pd(cloud->l1 + i, mul_pd(vdt, cloud->getVx1_pd(i))); // positionX tidbit
            store_pd(cloud->m1 + i, div_pd(mul_pd(vdt, load_pd(pFy)), vmass)); // velocityY tidbit
            store_pd(cloud->n1 + i, mul_pd(vdt, cloud->getVy1_pd(i))); // positionY tidbit

			// reset forces to zero:
			store_pd(pFx, set0_pd());
			store_pd(pFy, set0_pd());
		END_PARALLEL_FOR
        
		operate2(currentTime + dt/2.0);
		force2(currentTime + dt/2.0); // compute net force2
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k2 and l2 for entire cloud
            const doubleV vmass = load_pd(cloud->mass + i);

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;

            store_pd(cloud->k2 + i, div_pd(mul_pd(vdt, load_pd(pFx)), vmass)); // velocityX tidbit
            store_pd(cloud->l2 + i, mul_pd(vdt, cloud->getVx2_pd(i))); // positionX tidbit
            store_pd(cloud->m2 + i, div_pd(mul_pd(vdt, load_pd(pFy)), vmass)); // velocityY tidbit
            store_pd(cloud->n2 + i, mul_pd(vdt, cloud->getVy2_pd(i))); // positionY tidbit

			// reset forces to zero:
			store_pd(pFx, set0_pd());
			store_pd(pFy, set0_pd());
		END_PARALLEL_FOR

		operate3(currentTime + dt/2.0);
		force3(currentTime + dt/2.0); // compute net force3
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k3 and l3 for entire cloud
			const doubleV vmass = load_pd(cloud->mass + i);

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;

			store_pd(cloud->k3 + i, div_pd(mul_pd(vdt, load_pd(pFx)), vmass)); // velocityX tidbit
			store_pd(cloud->l3 + i, mul_pd(vdt, cloud->getVx3_pd(i))); // positionX tidbit
			store_pd(cloud->m3 + i, div_pd(mul_pd(vdt, load_pd(pFy)), vmass)); // velocityY tidbit
			store_pd(cloud->n3 + i, mul_pd(vdt, cloud->getVy3_pd(i))); // positionY tidbit
			
			// reset forces to zero:
			store_pd(pFx, set0_pd());
			store_pd(pFy, set0_pd());
		END_PARALLEL_FOR
        
		operate4(currentTime + dt);
		force4(currentTime + dt); // compute net force4
		BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static) // calculate k4 and l4 for entire cloud
			const doubleV vmass = load_pd(cloud->mass + i);

			// assign force pointers:
			double * const pFx = cloud->forceX + i;
			double * const pFy = cloud->forceY + i;
            
			store_pd(cloud->k4 + i, div_pd(mul_pd(vdt, load_pd(pFx)), vmass)); // velocityX tidbit
			store_pd(cloud->l4 + i, mul_pd(vdt, cloud->getVx4_pd(i))); // positionX tidbit
			store_pd(cloud->m4 + i, div_pd(mul_pd(vdt, load_pd(pFy)), vmass)); // velocityY tidbit
			store_pd(cloud->n4 + i, mul_pd(vdt, cloud->getVy4_pd(i))); // positionY tidbit
			
			// reset forces to zero:
			store_pd(pFx, set0_pd());
			store_pd(pFy, set0_pd());
		END_PARALLEL_FOR

        // Calculate next position and next velocity for entire cloud.
        BEGIN_PARALLEL_FOR(i, e, numParticles, DOUBLE_STRIDE, static)
			const doubleV vk1 = load_pd(cloud->k1 + i);
			const doubleV vk2 = load_pd(cloud->k2 + i);
			const doubleV vk3 = load_pd(cloud->k3 + i);
			const doubleV vk4 = load_pd(cloud->k4 + i);

			const doubleV vl1 = load_pd(cloud->l1 + i);
			const doubleV vl2 = load_pd(cloud->l2 + i);
			const doubleV vl3 = load_pd(cloud->l3 + i);
			const doubleV vl4 = load_pd(cloud->l4 + i);

			const doubleV vm1 = load_pd(cloud->m1 + i);
			const doubleV vm2 = load_pd(cloud->m2 + i);
			const doubleV vm3 = load_pd(cloud->m3 + i);
			const doubleV vm4 = load_pd(cloud->m4 + i);

			const doubleV vn1 = load_pd(cloud->n1 + i);
			const doubleV vn2 = load_pd(cloud->n2 + i);
			const doubleV vn3 = load_pd(cloud->n3 + i);
			const doubleV vn4 = load_pd(cloud->n4 + i);

			// calculate next positions and velocities:
			plusEqual_pd(cloud->Vx + i, da(vk1, vk2, vk3, vk4));
			plusEqual_pd(cloud->x + i, da(vl1, vl2, vl3, vl4));
			plusEqual_pd(cloud->Vy + i, da(vm1, vm2, vm3, vm4));
			plusEqual_pd(cloud->y + i, da(vn1, vn2, vn3, vn4));
		END_PARALLEL_FOR

		currentTime += dt;
	}
}

// 4th order Runge-Kutta full timstep parameter change.
inline const doubleV Runge_Kutta4::da(const doubleV a1, const doubleV a2, 
                                      const doubleV a3, const doubleV a4) {
    // (a1 + 2.0*(a2 + a3) + a4)/6.0
    return div_pd(add_pd(a1, add_pd(mul_pd(add_pd(a2, a3), 2.0), a4)), 6.0);
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
