/*===- ThermalForceLocalized.h - libSimulation -================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef THERMALFORCELOCALIZED_H
#define THERMALFORCELOCALIZED_H

#include "Force.h"
#include "mtrand.h"	//MT header
#include "VectorCompatibility.h"

class ThermalForceLocalized : public Force
{	
public:
	MTRand mt;
	ThermalForceLocalized(Cloud *myCloud, double thermRed1, double thermRed2, double specifiedRadius);	//overloaded constructor
	~ThermalForceLocalized() {} //destructor

//public variables:
	double heatVal1;
	double heatVal2;

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

    void writeForce(fitsfile *file, int *error);
    void readForce(fitsfile *file, int *error);

private:
//private variables:
	double heatingRadius;

//private functions:
	void force(const unsigned int currentParticle, const double displacementX, const double displacementY);
	void force(const unsigned int currentParticle, const __m128d displacementX, const __m128d displacementY);
};

#endif /* THERMALFORCELOCALIZED_H */
