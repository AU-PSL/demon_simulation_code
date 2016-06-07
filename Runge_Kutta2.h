/**
* @file  Runge_Kutta2.h
* @brief Defines the data and methods of the Runge_Kutta2 class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef RUNGA_KUTTA2_H
#define RUNGA_KUTTA2_H

#include "Integrator.h"

class Runge_Kutta2 : public Integrator {
public:
	Runge_Kutta2(Cloud * const C, const ForceArray &FA, 
                 const double timeStep, const double startTime);
	~Runge_Kutta2() {}

	virtual void moveParticles(const double endTime);
    
protected:
	void operate1(const double currentTime) const; // rk substep 1
	void operate2(const double currentTime) const; // rk substep 2
    
	void force1(const double currentTime) const; // rk substep 1
	void force2(const double currentTime) const; // rk substep 2
};

#endif // RUNGE_KUTTA2_H
