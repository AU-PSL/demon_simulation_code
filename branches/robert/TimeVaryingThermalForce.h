/*===- TimeVaryingThermalForce.h - libSimulation -==============================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef TIMEVARYINGTHERMALFORCE_H
#define TIMEVARYINGTHERMALFORCE_H

#include "ThermalForce.h"

class TimeVaryingThermalForce1D : public ThermalForce1D
{
public:
	TimeVaryingThermalForce1D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingThermalForce1D() {} //destructor

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
//private variables
	double heatValScale;
	double heatValOffset;

//private functions:
	const double calculateHeatVal(const double currentTime) const;
};

class TimeVaryingThermalForce2D : public ThermalForce2D TimeVaryingThermalForce1D
{
public:
	TimeVaryingThermalForce2D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingThermalForce2D() {} //destructor
};

class TimeVaryingThermalForce3D : public ThermalForce3D TimeVaryingThermalForce2D
{
public:
	TimeVaryingThermalForce3D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingThermalForce3D() {} //destructor
};

#endif /* TIMEVARYINGTHERMALFORCE_H */
