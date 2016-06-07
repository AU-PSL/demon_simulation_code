/**
* @file  Integrator.cpp
* @class Integrator Integrator.h
*
* @brief Provides support for integration methods
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#include "Integrator.h"
#include "CacheOperator.h"
#include <cmath>
#include <limits>

#ifdef DISPATCH_QUEUES
#define BLOCK_VALUE_TIME currTimeStep    //!< ??UNKNOWN??
#define BLOCK_VALUE_DIST currDist        //!< ??UNKNOWN??
#else
#define BLOCK_VALUE_TIME currentTimeStep //!< ??UNKNOWN??
#define BLOCK_VALUE_DIST currentDist     //!< ??UNKNOWN??
#endif

/**
* @brief Constructor for the Integrator class
*
* @param[in] C         Cloud object
* @param[in] FA        Array of forces
* @param[in] timeStep  Simulation time step
* @param[in] startTime Simulation start time
**/
Integrator::Integrator(Cloud * const C, const ForceArray &FA,
                       const double timeStep, double startTime)
: currentTime(startTime), cloud(C), forces(FA), init_dt(timeStep),
operations({{new CacheOperator(C)}})
SEMAPHORES_MALLOC(1) {
    SEMAPHORES_INIT(1);
}

/**
* @brief Destructor for the Integrator class
**/
Integrator::~Integrator() {
	for (Operator * const opt : operations)
		delete opt;
    
    SEMAPHORES_FREE(1);
}

/**
* @brief Reduces timestep if partices within a given distance.
*
* @details If particle spacing is less than the specified distance reduce timestep by a
* 		   factor of 10 and recheck with disance reduced by a factor of 10. Once all
*          particle spacings are outside the specified distance use the current 
*          timestep. This allows fine grain control of reduced timesteps
*
* @param[in] currentDist     The current distance..?
* @param[in] currentTImeStep The current simulation timestep
*
* @return The new timestep
*
* @bug When changing this over to AVX to simplify, change force methods to
* 	   triangle and block forces. Triangle forces cover the outter loop force. Block
*	   covers the inner loop forces. AVX specific differences would go into that 
*	   section.
**/
const double Integrator::modifyTimeStep(float currentDist, double currentTimeStep) const {
	#ifdef __AVX__
	#error "Integrator::modifyTimeStep does not fully support AVX."
	#endif
		const cloud_index numParticles = cloud->n;
	    
	#ifdef DISPATCH_QUEUES
	    // You cannot block capture method arguments. Store these values in to non
	    // const block captured variables. The correct names subsequently used by
	    // BLOCK_VALUE_DIST and BLOCK_VALUE_TIME
	    __block float currDist = currentDist;
	    __block double currTimeStep = currentTimeStep;
	#endif
	    
		// Loop through entire cloud, or until reduction occurs. Reset innerIndex 
	    // after each loop iteration.
	#ifdef DISPATCH_QUEUES
		const cloud_index outerLoop = numParticles;
	#else
		const cloud_index outerLoop = numParticles - 1;
	#endif

    BEGIN_PARALLEL_FOR(outerIndex, e, outerLoop, FLOAT_STRIDE, dynamic)
		// calculate separation distance b/t adjacent elements:
		const floatV outPosX = loadFloatVector(cloud->x + outerIndex);
		const floatV outPosY = loadFloatVector(cloud->y + outerIndex);
	
        // seperation (a1 - a2, a3 - a4, a1 - a3, a2 - a4)
        floatV sepx = _mm_hsub_ps(outPosX, _mm_shuffle_ps(outPosX, outPosX, _MM_SHUFFLE(1, 2, 0, 3)));
        floatV sepy = _mm_hsub_ps(outPosY, _mm_shuffle_ps(outPosY, outPosY, _MM_SHUFFLE(1, 2, 0, 3)));
        tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
		
        // seperation (a1 - a2, a3 - a4, a1 - a4, a2 - a3) The lower half 
        // operation is a repeat of the lower half of the above.
        sepx = _mm_hsub_ps(outPosX, _mm_shuffle_ps(outPosX, outPosX, _MM_SHUFFLE(1, 3, 0, 2)));
        sepy = _mm_hsub_ps(outPosY, _mm_shuffle_ps(outPosY, outPosY, _MM_SHUFFLE(1, 3, 0, 2)));
		tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
        
		// Calculate separation distance b/t nonadjacent elements:
		for (cloud_index innerIndex = outerIndex + FLOAT_STRIDE; innerIndex < numParticles; innerIndex += FLOAT_STRIDE) {
			const floatV inPosX = loadFloatVector(cloud->x + innerIndex);
			const floatV inPosY = loadFloatVector(cloud->y + innerIndex);
			
            // seperation (a1 - b1, a2 - b2, a3 - b3, a4 - b4) 
			sepx = outPosX - inPosX;
			sepy = outPosY - inPosY;
			tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
			
            // seperation (a1 - b2, a2 - b3, a3 - b4, a4 - b5)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(3, 2, 1, 0));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(3, 2, 1, 0));
			tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
			
            // seperation (a1 - b3, a2 - b4, a3 - b1, a4 - b2)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(0, 3, 2, 1));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(0, 3, 2, 1));
			tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
			
            // seperation (a1 - b4, a2 - b1, a3 - b2, a4 - b3)
			sepx = outPosX - _mm_shuffle_ps(inPosX, inPosX, _MM_SHUFFLE(1, 0, 3, 2));
			sepy = outPosY - _mm_shuffle_ps(inPosY, inPosY, _MM_SHUFFLE(1, 0, 3, 2));
			tryToReduceTimeStep(sepx, sepy, BLOCK_VALUE_DIST, BLOCK_VALUE_TIME);
		}
	END_PARALLEL_FOR

    return BLOCK_VALUE_TIME;
}

/**
* @brief Reduces timestep if particles are within distance
*
* @param[in]     sepX      Particle separation in x-direction
* @param[in]     sepY 	   Particle separation in y-direction
* @param[in,out] distance  Distance to check if less than
* @param[in,out] time      Current simulation timestep
**/
inline void Integrator::tryToReduceTimeStep(const floatV sepx, const floatV sepy, float &distance, double &time) const {
	// If particles are too close, reduce time step:
	while (isWithInDistance(sepx, sepy, distance)) {
		// Only one thread should modify the distance and timesStep at a time.
		SEMAPHORE_WAIT(0)
		if (isWithInDistance(sepx, sepy, distance)) {
			distance /= 10.0f;
			time /= 10.0f;
		}
		SEMAPHORE_SIGNAL(0)
	}
}

/**
* @brief Loads a float vector with particle locations
*
* @param[in] x Particle location
**/
inline floatV Integrator::loadFloatVector(double * const x) {
#ifdef __AVX__
    return _mm256_set_ps((float)x[0], (float)x[1], (float)x[2], (float)x[3],
                         (float)x[4], (float)x[5], (float)x[6], (float)x[7]);
#else
    return _mm_set_ps((float)x[0], (float)x[1], (float)x[2], (float)x[3]);
#endif
}

/**
* @brief Checks if any particles are within a certain distance
*
* @param[in] a    Vector of particle separations in x-direction
* @param[in] b    Vector of particle separations in y-direction
* @param[in] dist Distance to check
* 
* @return True if particle is within the set distance
**/
inline bool Integrator::isWithInDistance(const floatV a, const floatV b, const float dist) {
	return (bool)movemask_ps(cmple_ps(sqrt_ps(add_ps(mul_ps(a, a), mul_ps(b, b))), set1_ps(dist)));
}
