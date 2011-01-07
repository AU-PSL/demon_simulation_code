/*===- RotationalForce.h - libSimulation -======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef ROTATIONALFORCE_H
#define ROTATIONALFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"
#include <iostream> //necessary for cout statement in force1_1D

class RotationalForce : public Force
{	
public:
	RotationalForce(Cloud * const myCloud, const double rmin, const double rmax, const double rotConst);
	~RotationalForce() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1_1D(const double currentTime); //rk substep 1
	void force2_1D(const double currentTime); //rk substep 2
	void force3_1D(const double currentTime); //rk substep 3
	void force4_1D(const double currentTime); //rk substep 4
	//RotationalForce currently undefined for 1D.
	// Calling above functions merely exits program.

	void force1_2D(const double currentTime); 
	void force2_2D(const double currentTime); 
	void force3_2D(const double currentTime); 
	void force4_2D(const double currentTime); 

	void force1_3D(const double currentTime); 
	void force2_3D(const double currentTime); 
	void force3_3D(const double currentTime); 
	void force4_3D(const double currentTime); 

	void writeForce(fitsfile *file, int *error, const int dimension) const;
	void readForce(fitsfile *file, int *error, const int dimension);

private:
//private variables:
	double innerRad;
	double outerRad;
	double rotationalConst;

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
};

#endif /* ROTATIONALFORCE_H */
