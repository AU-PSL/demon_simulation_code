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

typedef unsigned char force_index;
typedef long force_flags;

// Binary assignments for the bit-packed FORCES keyword in Fits file:
enum ForceFlag {
	ConfinementForceFlag = 1,          //00000000001
	DragForceFlag = 2,                 //00000000010
	ShieldedCoulombForceFlag = 4,      //00000000100
	RectConfinementForceFlag = 8,      //00000001000
	ThermalForceFlag = 16,             //00000010000
	ThermalForceLocalizedFlag = 32,    //00000100000
	DrivingForceFlag = 64,             //00001000000
	RotationalForceFlag = 128,         //00010000000
	TimeVaryingDragForceFlag = 256,    //00100000000
	TimeVaryingThermalForceFlag = 512, //01000000000
	MagneticForceFlag = 1024           //10000000000
};

class Force
{
public:
	Cloud * const cloud;

	Force(Cloud * const myCloud) : cloud(myCloud) {} 
	virtual ~Force() {} //implementation of virtual destructor

	//Note: currentTime parameter necessary for DrivingForce, unused in others
	virtual void force1(const double currentTime)=0; //rk substep 1
	virtual void force2(const double currentTime)=0; //rk substep 2
	virtual void force3(const double currentTime)=0; //rk substep 3
	virtual void force4(const double currentTime)=0; //rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const=0; //output force information to file:
	virtual void readForce(fitsfile * const file, int * const error)=0;        //read force information from file:
};

#endif /* FORCE_H */
