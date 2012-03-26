/*===- Integrator.cpp - libSimulation -=========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#include "Integrator.h"
#include "CacheOperator.h"
#include <cmath>
#include <limits>

Integrator::Integrator(Cloud * const C, const ForceArray &FA,
                       const double timeStep, double startTime)
: currentTime(startTime), cloud(C), forces(FA), init_dt(timeStep),
operations({{new CacheOperator(C)}})
SEMAPHORES_MALLOC(1) {
    SEMAPHORES_INIT(1);
}

Integrator::~Integrator() {
	for (Operator * const opt : operations)
		delete opt;
    
    SEMAPHORES_FREE(1);
}

// If particle spacing is less than the specified distance reduce timestep by a
// factor of 10 and recheck with disance reduced by a factor of 10. Once all
// particle spacings are outside the specified distance use the current 
// timestep. This allows fine grain control of reduced timesteps.
const double Integrator::modifyTimeStep(float currentDist, double currentTimeStep) const {
#ifdef __AVX__
#error "Integartor::modifyTimeStep does not fully support AVX."
#endif
	// set constants:	
	const cloud_index numPar = cloud->n;
	const float redFactor = 10.0f;
    
#ifdef DISPATCH_QUEUES
    // You cannot block capture method arguments. Store these values in to non
    // const block captured variables. The correct names supsequently used by
    // BLOCK_VALUE_DIST and BLOCK_VALUE_TIME
    __block float currDist = currentDist;
    __block double currTimeStep = currentTimeStep;
#endif
    
	// Loop through entire cloud, or until reduction occures. Reset innerIndex 
    // after each loop iteration.
#ifdef DISPATCH_QUEUES
	const cloud_index outerLoop = numPar;
#else
	const cloud_index outerLoop = numPar - 1;
#endif
    BEGIN_PARALLEL_FOR(outerIndex, e, outerLoop, FLOAT_STRIDE, dynamic)
		// caculate separation distance b/t adjacent elements:
		const floatV outPosX = loadFloatVector(cloud->x + outerIndex);
		const floatV outPosY = loadFloatVector(cloud->y + outerIndex);
	
        // seperation (a1 - a2, a3 - a4, a1 - a3, a2 - a4)
        floatV sepx = _mm_hsub_ps(outPosX, _mm_shuffle_ps(outPosX, outPosX, _MM_SHUFFLE(1, 2, 0, 3)));
        floatV sepy = _mm_hsub_ps(outPosY, _mm_shuffle_ps(outPosY, outPosY, _MM_SHUFFLE(1, 2, 0, 3)));
        
		// If particles are too close, reduce time step:
        while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
            // Only one thread should modify the distance and timesStep at a time.
            SEMAPHORE_WAIT(0)
            if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
                BLOCK_VALUE_DIST /= redFactor;
                BLOCK_VALUE_TIME /= redFactor;
            }
            SEMAPHORE_SIGNAL(0)
        }
		
        // seperation (a1 - a2, a3 - a4, a1 - a4, a2 - a3) The lower half 
        // operation is a repeat of the lower half of the above.
        sepx = _mm_hsub_ps(outPosX, _mm_shuffle_ps(outPosX, outPosX, _MM_SHUFFLE(1, 3, 0, 2)));
        sepy = _mm_hsub_ps(outPosY, _mm_shuffle_ps(outPosY, outPosY, _MM_SHUFFLE(1, 3, 0, 2)));
	
		// If particles are too close, reduce time step:
		while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
			// Only one thread should modify the distance and timesStep at a time.
			SEMAPHORE_WAIT(0)
			if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
				BLOCK_VALUE_DIST /= redFactor;
				BLOCK_VALUE_TIME /= redFactor;
			}
			SEMAPHORE_SIGNAL(0)
		}
        
		// Calculate separation distance b/t nonadjacent elements:
		for (cloud_index innerIndex = outerIndex + FLOAT_STRIDE; innerIndex < numPar; innerIndex += FLOAT_STRIDE) {
			const floatV inPosX = loadFloatVector(cloud->x + innerIndex);
			const floatV inPosY = loadFloatVector(cloud->y + innerIndex);
			
            // seperation (a1 - b1, a2 - b2, a3 - b3, a4 - b4) 
			sepx = outPosX - inPosX;
			sepy = outPosY - inPosY;
			
			// If particles are too close, reduce time step:
			while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
				// Only one thread should modify the distance and timesStep at a time.
				SEMAPHORE_WAIT(0)
				if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
					BLOCK_VALUE_DIST /= redFactor;
					BLOCK_VALUE_TIME /= redFactor;
				}
				SEMAPHORE_SIGNAL(0)
			}
			
            // seperation (a1 - b2, a2 - b3, a3 - b4, a4 - b5)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(3, 2, 1, 0));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(3, 2, 1, 0));
			
			// If particles are too close, reduce time step:
			while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
				// Only one thread should modify the distance and timesStep at a time.
				SEMAPHORE_WAIT(0)
				if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
					BLOCK_VALUE_DIST /= redFactor;
					BLOCK_VALUE_TIME /= redFactor;
				}
				SEMAPHORE_SIGNAL(0)
			}
			
            // seperation (a1 - b3, a2 - b4, a3 - b1, a4 - b2)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(0, 3, 2, 1));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(0, 3, 2, 1));
			
			// If particles are too close, reduce time step:
			while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
				// Only one thread should modify the distance and timesStep at a time.
				SEMAPHORE_WAIT(0)
				if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
					BLOCK_VALUE_DIST /= redFactor;
					BLOCK_VALUE_TIME /= redFactor;
				}
				SEMAPHORE_SIGNAL(0)
			}
			
            // seperation (a1 - b4, a2 - b1, a3 - b2, a4 - b3)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(1, 0, 3, 2));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(1, 0, 3, 2));
			
			// If particles are too close, reduce time step:
			while (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
				// Only one thread should modify the distance and timesStep at a time.
				SEMAPHORE_WAIT(0)
				if (isWithInDistance(sepx, sepy, BLOCK_VALUE_DIST)) {
					BLOCK_VALUE_DIST /= redFactor;
					BLOCK_VALUE_TIME /= redFactor;
				}
				SEMAPHORE_SIGNAL(0)
			}
		}
	END_PARALLEL_FOR

    return BLOCK_VALUE_TIME;
}

inline floatV Integrator::loadFloatVector(double * const x) {
#ifdef __AVX__
    return _mm256_set_ps((float)x[0], (float)x[1], (float)x[2], (float)x[3],
                         (float)x[4], (float)x[5], (float)x[6], (float)x[7]);
#else
    return _mm_set_ps((float)x[0], (float)x[1], (float)x[2], (float)x[3]);
#endif
}

inline bool Integrator::isWithInDistance(const floatV a, const floatV b, const float dist) {
	return (bool)movemask_ps(cmple_ps(sqrt_ps(add_ps(mul_ps(a, a), mul_ps(b, b))), set1_ps(dist)));
}
