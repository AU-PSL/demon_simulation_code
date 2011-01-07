/*===- RectConfinementForce.h - libSimulation -=================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef RECTCONFINEMENTFORCE_H
#define RECTCONFINEMENTFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class RectConfinementForce : public Force
{
public:
	RectConfinementForce(Cloud * const myCloud, double confineConstX, double confineConstY, double confineConstZ); //confinement consts must be positive!
	~RectConfinementForce() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1_1D(const double currentTime); //rk substep 1
	void force2_1D(const double currentTime); //rk substep 2
	void force3_1D(const double currentTime); //rk substep 3
	void force4_1D(const double currentTime); //rk substep 4

	void force1_2D(const double currentTime); 
	void force2_2D(const double currentTime); 
	void force3_2D(const double currentTime); 
	void force4_2D(const double currentTime); 

	void force1_3D(const double currentTime); 
	void force2_3D(const double currentTime); 
	void force3_3D(const double currentTime); 
	void force4_3D(const double currentTime); 

	void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//private variables:
	double confineX;
	double confineY;
	double confineZ;

//private functions:
	void force1D(const unsigned int currentParticle, const __m128d currentPositionX);
	void force2D(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY);
	void force3D(const unsigned int currentParticle, const __m128d currentPositionX, const __m128d currentPositionY, const __m128d currentPositionZ);
};

#endif /* RECTCONFINEMENTFORCE_H */
