/*===- ThermalForce.h - libSimulation -=========================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef THERMALFORCE_H
#define THERMALFORCE_H

#include "Force.h"
#include "mtrand.h"	//MT header

class ThermalForce : public Force
{	
public:
	MTRand mt;
	ThermalForce(Cloud *myCloud, double redFactor);	//overloaded constructor
	~ThermalForce() {} //destructor

//public variables:
	double heatVal;

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); //rk substep 1
	virtual void force2(const double currentTime); //rk substep 2
	virtual void force3(const double currentTime); //rk substep 3
	virtual void force4(const double currentTime); //rk substep 4

	virtual void writeForce(fitsfile *file, int *error);
	virtual void readForce(fitsfile *file, int *error);

private:
//private functions:
	void force(const unsigned int currentParticle);
};

#endif /* THERMALFORCE_H */
