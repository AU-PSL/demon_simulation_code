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
#include "VectorCompatibility.h"
#include <cmath>
#include <limits>

Integrator::Integrator(Cloud * const myCloud, Force ** const forces, const force_index forcesSize,
                       const double timeStep, double startTime)
: cloud(myCloud), theForce(forces), numForces(forcesSize), init_dt(timeStep), currentTime(startTime), 
numOperators(1), operations(new Operator*[numOperators])
{
	// Operators are order dependent.
	operations[0] = new CacheOperator(cloud);
}

Integrator::~Integrator()
{
	for (operator_index i = 0; i < numOperators; i++)
		delete operations[i];
	delete[] operations;
}

/*------------------------------------------------------------------------------
 * If a particle spacing is less than the specified distance reduce timestep by a
 * factor of 10 and recheck with disance reduced by a factor of 10. Once all
 * particle spacings are outside the specified distance use the current timestep.
 * This allows fine grain control of reduced timesteps.
 ------------------------------------------------------------------------------*/
const double Integrator::modifyTimeStep(cloud_index outerIndex, cloud_index innerIndex, const double currentDist, 
                                        const double currentTimeStep) const
{
	// set constants:	
	const cloud_index numPar = cloud->n;
	const __m128d distv = _mm_set1_pd(currentDist);
	const double redFactor = 10.0;
    
	// loop through entire cloud, or until reduction occures. Reset innerIndex after each loop iteration.
	for (cloud_index e = numPar - 1; outerIndex < e; outerIndex += 2, innerIndex = outerIndex + 2)
	{
		// caculate separation distance b/t adjacent elements:
		const double sepx = cloud->x[outerIndex] - cloud->x[outerIndex + 1];
		const double sepy = cloud->y[outerIndex] - cloud->y[outerIndex + 1];
        
		// if particles too close, reduce time step:
		if (sqrt(sepx*sepx + sepy*sepy) <= currentDist)
			return modifyTimeStep(outerIndex, innerIndex, currentDist/redFactor, currentTimeStep/redFactor);
        
		// load positions into vectors:
		const __m128d vx1 = cloud->getx1_pd(outerIndex);	// x vector
		const __m128d vy1 = cloud->gety1_pd(outerIndex);	// y vector
        
		// calculate separation distance b/t nonadjacent elements:
		for (; innerIndex < numPar; innerIndex += 2)
		{
			// assign position pointers:
			const double * const px2 = cloud->x + innerIndex;
			const double * const py2 = cloud->y + innerIndex;
            
			// calculate j,i and j+1,i+1 separation distances:
			__m128d vx2 = vx1 - _mm_load_pd(px2);
			__m128d vy2 = vy1 - _mm_load_pd(py2);
            
			// check separation distances against dist:
			if (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), distv))	// if either are too close, reduce time step
				return modifyTimeStep(outerIndex, innerIndex, currentDist/redFactor, currentTimeStep/redFactor);
            
			// calculate j,i+1 and j+1,i separation distances:
			vx2 = vx1 - _mm_loadr_pd(px2);
			vy2 = vy1 - _mm_loadr_pd(py2);
            
			// check separation distances against dist:
			if (isLessThanOrEqualTo(_mm_sqrt_pd(vx2*vx2 + vy2*vy2), distv))	// if either are too close, reduce time step
				return modifyTimeStep(outerIndex, innerIndex, currentDist/redFactor, currentTimeStep/redFactor);
		}
	}
    
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
