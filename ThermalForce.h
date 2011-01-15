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
#include "mtrand.h" //MT header

class ThermalForce1D : public Force
{
public:
	ThermalForce1D(Cloud * const myCloud, const double redFactor);
	~ThermalForce1D() {}; //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); //rk substep 1
	virtual void force2(const double currentTime); //rk substep 2
	virtual void force3(const double currentTime); //rk substep 3
	virtual void force4(const double currentTime); //rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	virtual void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//private variables:
	MTRand mt;

//private functions:
	void force(const unsigned int currentParticle);

protected:
//protected variables:
	double heatVal;
};

class ThermalForce2D : public ThermalForce1D
{
public:
	ThermalForce2D(Cloud * const myCloud, const doudble redFactor);
	~ThermalForce2D() {};

private:
//private functions:
	void force(const unsigned int currentParticle);
};

class ThermalForce3D : public ThermalForce1D
{
public:
	ThermalForce3D(Cloud * const myCloud, const doudble redFactor);
	~ThermalForce3D() {};

private:
//private functions:
	void force(const unsigned int currentParticle);
};

#endif /* THERMALFORCE_H */
