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
	ThermalForceLocalized(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius);	//overloaded constructor
	~ThermalForceLocalized() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

    void writeForce(fitsfile * const file, int * const error);
    void readForce(fitsfile * const file, int * const error);

private:
//private variables:
    MTRand mt;
	double heatingRadius;
    double heatVal1;
	double heatVal2;

//private functions:
	void force(const unsigned int currentParticle, const double displacementX, const double displacementY);
	void force(const unsigned int currentParticle, const __m128d displacementX, const __m128d displacementY);
};

#endif /* THERMALFORCELOCALIZED_H */
