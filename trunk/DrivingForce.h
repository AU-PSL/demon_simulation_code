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

class DrivingForce : public Force {
public:
	DrivingForce(Cloud * const C, const double dampConst, const double amp, const double drivingShift)
	: Force(C), amplitude(amp), driveConst(-dampConst), shift(drivingShift) {}
	~DrivingForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double amplitude, driveConst, shift; // [N], [m^2], [m]
	static const double waveNum; // [m^-1]
	static const double angFreq; // [rad*Hz]

	void force(const cloud_index currentParticle, const doubleV currentTime, const doubleV currentPositionX);
};

#endif // DRIVINGFORCE_H
