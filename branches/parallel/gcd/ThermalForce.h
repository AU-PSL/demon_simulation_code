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
#include "mtrand.h"	// MT header

class ThermalForce : public Force
{	
public:
	ThermalForce(Cloud * const myCloud, const double redFactor);
	~ThermalForce(); // destructor

// public functions:
	// Note: currentTime parameter is necessary (due to parent class) but unused
	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
// private variables:
	MTRand mt;
	dispatch_semaphore_t semaphore;

// private functions:
	void force(const cloud_index currentParticle);

protected:
// protected variables:
	double heatVal;
};

#endif // THERMALFORCE_H
