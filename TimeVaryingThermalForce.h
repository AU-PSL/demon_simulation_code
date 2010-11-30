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

class TimeVaryingThermalForce : public ThermalForce
{	
public:
	TimeVaryingThermalForce(Cloud *myCloud, double scale, double offset);	//overloaded constructor
	~TimeVaryingThermalForce() {} //destructor

//public variables:
	double heatValScale;
    double heatValOffset;

//public functions:
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

    void writeForce(fitsfile *file, int *error);
    void readForce(fitsfile *file, int *error);
};

#endif /* TIMEVARYINGTHERMALFORCE_H */
