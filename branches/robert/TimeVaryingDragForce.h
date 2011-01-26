/*===- TimeVaryingDragForce.h - libSimulation -=================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef TIMEVARYINGDRAGFORCE_H
#define TIMEVARYINGDRAGFORCE_H

#include "DragForce.h"

class TimeVaryingDragForce1D : public DragForce1D
{
public:
	TimeVaryingDragForce1D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingDragForce1D() {} //destructor

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

protected:
//protected functions:
	const double calculateGamma(const double currentTime) const;

private:
//private variables:
	double scaleConst;  //[s^-2]
	double offsetConst; //[s^-1]
};

class TimeVaryingDragForce2D : public TimeVaryingDragForce1D, public DragForce2D
{
public:
	TimeVaryingDragForce2D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingDragForce2D() {} //destructor

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4
};

class TimeVaryingDragForce3D : public TimeVaryingDragForce2D, public DragForce3D
{
public:
	TimeVaryingDragForce3D(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingDragForce3D() {} //destructor

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4
};

#endif /* TIMEVARYINGDRAGFORCE_H */
