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
#include <array>

class Integrator {
public:
    Integrator(Cloud * const myCloud, Force ** const forces, const force_index forcesSize,
               const double timeStep, double startTime);
    virtual ~Integrator();
    
    // public variables:
	Cloud * const cloud; // pointer to cloud object
	Force ** const theForce; // pointer to Force object
	const force_index numForces;
	const double init_dt; // store initial time step
	double currentTime;
    
    virtual void moveParticles(const double endTime)=0;
    
protected:
    const std::array<Operator *, 1> operations;
    SEMAPHORES
    
    const double modifyTimeStep(float currentDist, double currentTimeStep) const;
	static __m128 loadFloatVector(double * const x);
	static bool isWithInDistance(const __m128 a, const __m128 b, const float dist);
};

#endif // INTEGRATOR_H
