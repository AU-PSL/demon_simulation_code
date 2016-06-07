/**
* @file  RectConfinementForce.h
* @brief Defines the data and methods of the RectConfinementForce class
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef RECTCONFINEMENTFORCE_H
#define RECTCONFINEMENTFORCE_H

#include "Force.h"

class RectConfinementForce : public Force {
public:

	/**
	* @brief Constructor for the RectConfinementForce class
	**/
	RectConfinementForce(Cloud * const C, double confineConstX, double confineConstY)
	: Force(C), confineX(confineConstX), confineY(confineConstY) {}
	
	/**
	* @brief Destructor for the RectConfinementForce class
	**/
	~RectConfinementForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4

	void writeForce(fitsfile * const file, int * const error) const;
	void readForce(fitsfile * const file, int * const error);

private:
	double confineX; //<! Strength of the confinement force in the x-direction
	double confineY; //<! Strength of the confinement force in the y-direction

	void force(const cloud_index currentParticle, const doubleV currentPositionX, const doubleV currentPositionY);
};

#endif // RECTCONFINEMENTFORCE_H
