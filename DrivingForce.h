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
	DrivingForce(Cloud *myCloud, double dampConst, double amp, double drivingShift);	//overloaded constructor
	~DrivingForce() {} //destructor

//public variables:
	double amplitude;	//[m]

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

    void writeForce(fitsfile *file, int *error);
    void readForce(fitsfile *file, int *error);

private:
//private variables:
	double driveConst;	//[m]
	double shift;		//[m]
	static const double waveNum;		//[1/m]
	static const double angFreq;		//[rad*Hz]

//private functions:
	void force(const unsigned int currentParticle, const __m128d currentTime, const __m128d currentPositionX);
};

#endif /* DRIVINGFORCE_H */
