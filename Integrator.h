/**
* @file  Integrator.h
* @brief Defines the data and methods of the Integrator class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Cloud.h"
#include "Force.h"
#include "Operator.h"
#include <array>

class Integrator {
public:
    Integrator(Cloud * const C, const ForceArray &FA,
               const double timeStep, double startTime);
    virtual ~Integrator();

	double currentTime;
    
    virtual void moveParticles(const double endTime)=0;
    
protected:
	Cloud * const cloud; // pointer to cloud object
	const ForceArray &forces;
	const double init_dt; // store initial time step
    const std::array<Operator * const, 1> operations;
    SEMAPHORES
    
    const double modifyTimeStep(float currentDist, double currentTimeStep) const;
	void tryToReduceTimeStep(const floatV sepx, const floatV sepy, float &distance, double &time) const;
	static floatV loadFloatVector(double * const x);
	static bool isWithInDistance(const floatV a, const floatV b, const float dist);
};

#endif // INTEGRATOR_H
