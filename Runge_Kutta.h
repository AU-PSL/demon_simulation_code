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
	Runge_Kutta(Cloud *myCloud, Force **forces, double timeStep, unsigned int forcesSize, double startTime);	//overloaded constructor

//public variables:
	Cloud *cloud;		//pointer to cloud object
	Force **theForce;	//pointer to Force object
	unsigned int numForces;
	double dt;		//current time step
	double init_dt;		//store initial time step
	double red_dt;		//store reduced time step
	double currentTime;

//public functions:
	//Input: double endTime
	//Preconditions: endTime > 0
	//Postconditions: Runge-Kutta algorithm complete; position, velocity, time updated.
	void moveParticles(double endTime);

private:
//private functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4
	void modifyTimeStep();
};

#endif /* RUNGE_KUTTA_H */
