/**
* @file  VertElectricForce.h
* @brief Defines the data and methods of the VertElectricForce class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef VERTELECTRICFORCE_H
#define VERTELECTRICFORCE_H

#include "Force.h"

class VertElectricForce : public Force {
public:
	VertElectricForce(Cloud * const C, double vertElectricField, double vertDecay)
	: Force(C), vertElectric(vertElectricField), vertDec(vertDecay) {}
	~VertElectricForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
	double vertElectric; //<! Strength of vertical electric force [V/m^2]
	double vertDec;      //<! Strength of the vertical decay factor [m]

	void force(const cloud_index currentParticle, const doubleV currentPositionY);
};

#endif // VERTELECTRICFORCE_H
