/*===- Runge_Kutta4.h - libSimulation -=========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef RUNGA_KUTTA4_H
#define RUNGA_KUTTA4_H

#include "Runge_Kutta2.h"
#include "Force.h"
#include "Operator.h"

class Runge_Kutta4 : public Runge_Kutta2 {
public:
	Runge_Kutta4(Cloud * const C, const ForceArray &FA, 
                 const double timeStep, const double startTime);
	~Runge_Kutta4() {}

// public functions:
	// Input: double endTime
	// Preconditions: endTime > 0
	// Postconditions: Runge-Kutta algorithm complete; position, velocity, time updated.
	void moveParticles(const double endTime);

private:
// private functions:
	void operate3(const double currentTime) const; // rk substep 3
	void operate4(const double currentTime) const; // rk substep 4
  
	void force3(const double currentTime) const; // rk substep 3
	void force4(const double currentTime) const; // rk substep 4
};

#endif // RUNGE_KUTTA4_H
