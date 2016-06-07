/**
* @file  RotationalForce.h
* @brief Defines the data and methods of the RotationalForce class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef ROTATIONALFORCE_H
#define ROTATIONALFORCE_H

#include "Force.h"

class RotationalForce : public Force {
public:

	/**
	* @brief Constructor for the RotationalForce class
	**/
	RotationalForce(Cloud * const C, const double rmin, const double rmax, const double rotConst)
	: Force(C), innerRad(rmin), outerRad(rmax), rotationalConst(rotConst) {}

	/**
	* @brief Destructor for the RotationalForce class
	**/
	~RotationalForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

private:
	double innerRad;		//<! F=0 if radius is less than this value [m]
	double outerRad;		//<! F=0 if radius is less than this value [m]
	double rotationalConst; //<! Strength of rotational force [N]

	void force(const cloud_index currentParticle, 
               const doubleV currentPositionX, 
               const doubleV currentPositionY);
};

#endif // ROTATIONALFORCE_H
