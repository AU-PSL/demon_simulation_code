/*===- Runge_Kutta.h - libSimulation -==========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef RUNGA_KUTTA_H
#define RUNGA_KUTTA_H

#include "Force.h"

class Operator;
typedef unsigned char operator_index;

class Runge_Kutta
{
public:
	Runge_Kutta(Cloud * const myCloud, Force ** const forces, const double timeStep, const force_index forcesSize, double startTime);
	~Runge_Kutta();

// public variables:
	Cloud * const cloud; // pointer to cloud object
	Force ** const theForce; // pointer to Force object
	const force_index numForces;
	const double init_dt; // store initial time step
	double currentTime;

// public functions:
	// Input: double endTime
	// Preconditions: endTime > 0
	// Postconditions: Runge-Kutta algorithm complete; position, velocity, time updated.
	void moveParticles(const double endTime);

private:
// private variables:
	const operator_index numOperators;
	Operator ** const operations;
	dispatch_queue_t queue;
	dispatch_semaphore_t sema;
    
// private functions:
	void operate1(const double currentTime) const; // rk substep 1
	void operate2(const double currentTime) const; // rk substep 2
	void operate3(const double currentTime) const; // rk substep 3
	void operate4(const double currentTime) const; // rk substep 4
  
	void force1(const double currentTime) const; // rk substep 1
	void force2(const double currentTime) const; // rk substep 2
	void force3(const double currentTime) const; // rk substep 3
	void force4(const double currentTime) const; // rk substep 4
	
	const double modifyTimeStep(const double currentDist, const double currentTimeStep) const;
	static bool isLessThanOrEqualTo(const __m128d a, const __m128d b);
};

#endif // RUNGE_KUTTA_H
