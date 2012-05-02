// This is a test
/*===- Force.h - libSimulation -================================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef FORCE_H
#define FORCE_H

#include "Cloud.h"
#include <vector>

class Force {
public:
	Cloud * const cloud;
	
	Force(Cloud * const C) : cloud(C) {} 
	virtual ~Force() {}

	virtual void force1(const double currentTime)=0; // rk substep 1
	virtual void force2(const double currentTime)=0; // rk substep 2
	virtual void force3(const double currentTime)=0; // rk substep 3
	virtual void force4(const double currentTime)=0; // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const=0;	// output force information to file
	virtual void readForce(fitsfile * const file, int * const error)=0;	// read force information from file
};
	
typedef std::vector<Force *> ForceArray;
typedef long force_flags;

// Binary assignments for the bit-packed FORCES keyword in Fits file:
enum ForceFlag : force_flags {
	ConfinementForceFlag = 1,          // 000000000001
	DragForceFlag = 2,                 // 000000000010
	ShieldedCoulombForceFlag = 4,      // 000000000100
	RectConfinementForceFlag = 8,      // 000000001000
	ThermalForceFlag = 16,             // 000000010000
	ThermalForceLocalizedFlag = 32,    // 000000100000
	DrivingForceFlag = 64,             // 000001000000
	RotationalForceFlag = 128,         // 000010000000
	TimeVaryingDragForceFlag = 256,    // 000100000000
	TimeVaryingThermalForceFlag = 512, // 001000000000
	MagneticForceFlag = 1024,          // 010000000000
	ConfinementForceVoidFlag = 2048    // 100000000000
};

#endif // FORCE_H
