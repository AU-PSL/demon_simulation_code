/*===- DrivingForce.h - libSimulation -=========================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef DRIVINGFORCE_H
#define DRIVINGFORCE_H

#include "Force.h"
#include "VectorCompatibility.h"

class DrivingForce : public Force
{
public:
	DrivingForce(Cloud * const myCloud, const double dampConst, const double amp, const double drivingShift); //overloaded constructor
	~DrivingForce() {} //destructor

//public functions:
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
	double amplitude;  //[N]
	double driveConst; //[m]
	double shift;      //[m]
	static const double waveNum; //[1/m]
	static const double angFreq; //[rad*Hz]

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentTime, const __m128d currentPositionX);
};

#endif /* DRIVINGFORCE_H */
