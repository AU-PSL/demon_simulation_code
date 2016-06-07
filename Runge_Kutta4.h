/**
* @file  Runge_Kutta4.h
* @brief Defines the data and methods of the Runge_Kutta4 class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef RUNGA_KUTTA4_H
#define RUNGA_KUTTA4_H

#include "Runge_Kutta2.h"

class Runge_Kutta4 : public Runge_Kutta2 {
public:
	Runge_Kutta4(Cloud * const C, const ForceArray &FA, 
                 const double timeStep, const double startTime);
	~Runge_Kutta4() {}

	void moveParticles(const double endTime);

private:
	void operate3(const double currentTime) const; // rk substep 3
	void operate4(const double currentTime) const; // rk substep 4
  
	void force3(const double currentTime) const; // rk substep 3
	void force4(const double currentTime) const; // rk substep 4
    
    static const doubleV da(const doubleV a1, const doubleV a2, 
                            const doubleV a3, const doubleV a4);
};

#endif // RUNGE_KUTTA4_H
