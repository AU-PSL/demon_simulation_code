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
#include "mtrand.h" //MT header
#include "VectorCompatibility.h"
#include <cmath>
#include <ctime>

class ThermalForceLocalized1D : public Force
{
public:
	ThermalForceLocalized1D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius);
	~ThermalForceLocalized1D() {}

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	void force1(const double currentTime); //rk substep 1
	void force2(const double currentTime); //rk substep 2
	void force3(const double currentTime); //rk substep 3
	void force4(const double currentTime); //rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
//private variables:
	MTRand mt;
	double heatingRadius;
	double heatVal1;
	double heatVal2;

//private functions:
	void force1D(const unsigned int currentParticle, const __m128d displacementX);
};

class ThermalForceLocalized2D : public ThermalForceLocalized1D
{
public:
	ThermalForceLocalized2D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius);
	~ThermalForceLocalized2D() {}

private:
//private functions:
	void force2D(const unsigned int currentParticle, const __m128d displacementX, const __m128d displacementY);
};

class ThermalForceLocalized3D : public ThermalForceLocalized2D
{
public:
	ThermalForceLocalized3D(Cloud * const myCloud, const double thermRed1, const double thermRed2, const double specifiedRadius);
	~ThermalForceLocalized3D() {}

private:
//private functions:
	void force3D(const unsigned int currentParticle, const __m128d displacementX, const __m128d displacementY, const __m128d displacementZ);
};

#endif /* THERMALFORCELOCALIZED_H */
