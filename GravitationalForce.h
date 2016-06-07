/**
* @file  GravitationalForce.h
* @brief Defines the data and methods of the GravitationalForce class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef GRAVITATIONALFORCE_H
#define GRAVITATIONALFORCE_H

#include "Force.h"

class GravitationalForce : public Force {
public:
	GravitationalForce(Cloud * const C, const double gravitationalField)
	: Force(C), gravitational(gravitationalField) {}
	~GravitationalForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4
	
	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

protected:
	double gravitational; //!< The strength of the gravitational force [m/s^2]

private:
	void force(const cloud_index currentParticle);
};

#endif // GRAVITATIONALFORCE_H
