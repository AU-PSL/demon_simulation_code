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

class TimeVaryingDragForce : public DragForce
{
public:
	TimeVaryingDragForce(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingDragForce() {} //destructor

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

	void writeForce(fitsfile *file, int *error, const int dimension) const;
	void readForce(fitsfile *file, int *error, const int dimension);

private:
//private variables:
	double scaleConst;  //[s^-2]
	double offsetConst; //[s^-1]

//private methods:
	const double calculateGamma(const double currentTime) const;
};

#endif /* TIMEVARYINGDRAGFORCE_H */
