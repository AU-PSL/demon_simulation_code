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

class Runge_Kutta
{
public:
	Runge_Kutta(Cloud * const myCloud, Force ** const forces, const double timeStep, const unsigned int forcesSize, double startTime);	//overloaded constructor

//public variables:
	Cloud * const cloud;		//pointer to cloud object
	Force ** const theForce;	//pointer to Force object
	const unsigned int numForces;
	const double init_dt;		//store initial time step
	const double red_dt;        //store reduced time step
	double currentTime;

//public functions:
	//Input: double endTime
	//Preconditions: endTime > 0
	//Postconditions: Runge-Kutta algorithm complete; position, velocity, time updated.
	void moveParticles(const double endTime);

private:
//private functions:
	void force1(const double currentTime) const; //rk substep 1
	void force2(const double currentTime) const; //rk substep 2
	void force3(const double currentTime) const; //rk substep 3
	void force4(const double currentTime) const; //rk substep 4
	const double modifyTimeStep() const;
};

#endif /* RUNGE_KUTTA_H */