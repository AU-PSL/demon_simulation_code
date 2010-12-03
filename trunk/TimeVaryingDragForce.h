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
using namespace std;

class TimeVaryingDragForce : public DragForce
{	
public:
	TimeVaryingDragForce(Cloud * const myCloud, const double scale, const double offset);	//overloaded constructor
	~TimeVaryingDragForce() {} //destructor

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4
	
	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

private:
//private variables:
    double scaleConst;	//[s^-2]
    double offsetConst; //[s^-1]
	
//private methods:
	const double calculateGamma(const double currentTime) const;
};

#endif /* TIMEVARYINGDRAGFORCE_H */
