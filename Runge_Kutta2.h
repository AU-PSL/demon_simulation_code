/*===- Runge_Kutta2.h - libSimulation -=========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef RUNGA_KUTTA2_H
#define RUNGA_KUTTA2_H

#include "Integrator.h"
#include "Force.h"
#include "Operator.h"

class Runge_Kutta2 : public Integrator {
public:
	Runge_Kutta2(Cloud * const myCloud, Force ** const forces, const force_index forcesSize, 
                 const double timeStep, const double startTime);
	~Runge_Kutta2() {}
    
    // public functions:
	// Input: double endTime
	// Preconditions: endTime > 0
	// Postconditions: Runge-Kutta algorithm complete; position, velocity, time updated.
	void moveParticles(const double endTime);
    
protected:
    // private functions:
	void operate1(const double currentTime) const; // rk substep 1
	void operate2(const double currentTime) const; // rk substep 2
    
	void force1(const double currentTime) const; // rk substep 1
	void force2(const double currentTime) const; // rk substep 2
};

#endif // RUNGE_KUTTA2_H
