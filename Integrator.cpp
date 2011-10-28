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
#include "Parallel.h"
#include <cmath>
#include <limits>

Integrator::Integrator(Cloud * const myCloud, Force ** const forces, const force_index forcesSize,
                       const double timeStep, double startTime)
: cloud(myCloud), theForce(forces), numForces(forcesSize), init_dt(timeStep), currentTime(startTime), 
numOperators(1), operations(new Operator*[numOperators]) {
	// Operators are order dependent.
	operations[0] = new CacheOperator(cloud);
}

Integrator::~Integrator() {
	begin_parallel_for(i, e, numOperators, 1)
		delete operations[i];
    end_parallel_for
	delete[] operations;
}

/*------------------------------------------------------------------------------
* If a particle spacing is less than the specified distance reduce timestep by a
* factor of 10 and recheck with disance reduced by a factor of 10. Once all
* particle spacings are outside the specified distance use the current timestep.
* This allows fine grain control of reduced timesteps.
------------------------------------------------------------------------------*/
const double Integrator::modifyTimeStep(double currentDist, double currentTimeStep) const {
	// set constants:	
	const cloud_index numPar = cloud->n;
	const double redFactor = 10.0;
    
	// loop through entire cloud, or until reduction occures. Reset innerIndex after each loop iteration.
	begin_parallel_for(outerIndex, e, cloud->n - 1, 2)
		// caculate separation distance b/t adjacent elements:
		const double sepx = cloud->x[outerIndex] - cloud->x[outerIndex + 1];
		const double sepy = cloud->y[outerIndex] - cloud->y[outerIndex + 1];
        
		// if particles too close, reduce time step:
		while (sqrt(sepx*sepx + sepy*sepy) <= currentDist)
            if (sqrt(sepx*sepx + sepy*sepy) <= currentDist) {
                currentDist /= redFactor;
                currentTimeStep /= redFactor;
            }
		
		// load positions into vectors:
		const __m128d vx1 = cloud->getx1_pd(outerIndex);	// x vector
		const __m128d vy1 = cloud->gety1_pd(outerIndex);	// y vector
        
		// calculate separation distance b/t nonadjacent elements:
		for (cloud_index innerIndex = outerIndex + 2; innerIndex < numPar; innerIndex += 2) {
			// assign position pointers:
			const double * const px2 = cloud->x + innerIndex;
			const double * const py2 = cloud->y + innerIndex;
            
			// calculate j,i and j+1,i+1 separation distances:
			__m128d vx2 = vx1 - _mm_load_pd(px2);
			__m128d vy2 = vy1 - _mm_load_pd(py2);
            
			// check separation distances against dist. If either are too close, reduce time step.
			while (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), _mm_set1_pd(currentDist)))
                // Retest condition to make sure a different thread hasn't already reduced.
                if (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), _mm_set1_pd(currentDist))) {
                    currentDist /= redFactor;
                    currentTimeStep /= redFactor;
                }
            
			// calculate j,i+1 and j+1,i separation distances:
			vx2 = vx1 - _mm_loadr_pd(px2);
			vy2 = vy1 - _mm_loadr_pd(py2);
            
			// check separation distances against dist. If either are too close, reduce time step.
			while (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), _mm_set1_pd(currentDist)))
                // Retest condition to make sure a different thread hasn't already reduced.
                if (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), _mm_set1_pd(currentDist))) {
                    currentDist /= redFactor;
                    currentTimeStep /= redFactor;
                }
		}
	end_parallel_for
    
	// reset time step:
	return currentTimeStep;
}

inline bool Integrator::isLessThanOrEqualTo(const __m128d a, const __m128d b) {
	__m128d comp = _mm_cmple_pd(a, b);
	
	double low, high;
	_mm_storel_pd(&low, comp);
	_mm_storeh_pd(&high, comp);
	
	return std::isnan(low) || std::isnan(high);
}
