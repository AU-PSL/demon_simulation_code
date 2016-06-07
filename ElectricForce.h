/**
* @file  ElectricForce.h
* @brief Defines the data and methods of the ElectricForce class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef ELECTRICFORCE_H
#define ELECTRICFORCE_H

#include "Force.h"

class ElectricForce : public Force {
public:
	ElectricForce(Cloud * const C, double electricField, double plasmaRad)
	: Force(C), electric(electricField), radius(plasmaRad) {}
	~ElectricForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4

	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

private:
	double electric; //!< Strength of the electric force [V/m^2]
	double radius; 	 //!< Decay constant of the electric force [m]

	void force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY);
};

#endif // ELECTRICFORCE_H
