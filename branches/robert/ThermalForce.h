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
	ThermalForce(Cloud * const myCloud, const double redFactor);	//overloaded constructor
	~ThermalForce() {} //destructor

//public functions:
	//Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1_1D(const double currentTime); //rk substep 1
	virtual void force2_1D(const double currentTime); //rk substep 2
	virtual void force3_1D(const double currentTime); //rk substep 3
	virtual void force4_1D(const double currentTime); //rk substep 4

	virtual void force1_2D(const double currentTime); 
	virtual void force2_2D(const double currentTime); 
	virtual void force3_2D(const double currentTime); 
	virtual void force4_2D(const double currentTime); 

	virtual void force1_3D(const double currentTime); 
	virtual void force2_3D(const double currentTime); 
	virtual void force3_3D(const double currentTime); 
	virtual void force4_3D(const double currentTime); 

	virtual void writeForce(fitsfile * const file, int * const error, const int dimension) const;
	virtual void readForce(fitsfile * const file, int * const error, const int dimension);

private:
//private variables:
	MTRand mt;

//private functions:
	void force1D(const unsigned int currentParticle);
	void force2D(const unsigned int currentParticle);
	void force3D(const unsigned int currentParticle);

protected:
//protected variables:
	double heatVal;
};

#endif /* THERMALFORCE_H */
