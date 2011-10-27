/*===- Integrator.h - libSimulation -===========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Cloud.h"
#include "Force.h"
#include "Operator.h"

class Integrator {
public:
    Integrator(Cloud * const myCloud, Force ** const forces, const force_index forcesSize,
               const double timeStep, double startTime);
    ~Integrator();
    
    // public variables:
	Cloud * const cloud; // pointer to cloud object
	Force ** const theForce; // pointer to Force object
	const force_index numForces;
	const double init_dt; // store initial time step
	double currentTime;
    
    virtual void moveParticles(const double endTime)=0;
    
protected:
    const operator_index numOperators;
	Operator ** const operations;
    
    const double modifyTimeStep(cloud_index outerIndex, cloud_index innerIndex, const double currentDist, 
                                const double currentTimeStep) const;
	static bool isLessThanOrEqualTo(const __m128d a, const __m128d b);
};

#endif // INTEGRATOR_H
