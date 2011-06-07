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

//TimeVaryingThermalForce1D inherits from ThermalForce3D (as opposed to ThermalForce1D) to circumvent
// the need for multiple inheritance in TimeVaryingThermalForce2D/3D.
class TimeVaryingThermalForce1D : public ThermalForce3D
{
public:
	TimeVaryingThermalForce1D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	virtual ~TimeVaryingThermalForce1D() {} //destructor

//public functions:
	virtual void force1(const double currentTime); //rk substep 1
	virtual void force2(const double currentTime); //rk substep 2
	virtual void force3(const double currentTime); //rk substep 3
	virtual void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

protected:
//protected functions:
	const double calculateHeatVal(const double currentTime) const;

private:
//private variables
	double heatValScale;
	double heatValOffset;
};

class TimeVaryingThermalForce2D : public TimeVaryingThermalForce1D
{
public:
	TimeVaryingThermalForce2D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	virtual ~TimeVaryingThermalForce2D() {} //destructor

//public functions:
	virtual void force1(const double currentTime);
	virtual void force2(const double currentTime);
	virtual void force3(const double currentTime);
	virtual void force4(const double currentTime);
};

class TimeVaryingThermalForce3D : public TimeVaryingThermalForce2D
{
public:
	TimeVaryingThermalForce3D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingThermalForce3D() {} //destructor

//public functions:
	void force1(const double currentTime);
	void force2(const double currentTime);
	void force3(const double currentTime);
	void force4(const double currentTime);
};

#endif /* TIMEVARYINGTHERMALFORCE_H */
