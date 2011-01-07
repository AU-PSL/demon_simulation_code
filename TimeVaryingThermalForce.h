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
	TimeVaryingThermalForce(Cloud * const myCloud, const double scale, const double offset); //overloaded constructor
	~TimeVaryingThermalForce() {} //destructor

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

	void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//private variables
	double heatValScale;
	double heatValOffset;

//private functions:
	const double calculateHeatVal(const double currentTime) const;
};

#endif /* TIMEVARYINGTHERMALFORCE_H */
