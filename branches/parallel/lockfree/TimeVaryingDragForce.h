/*===- TimeVaryingDragForce.h - libSimulation -=================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
* 
*===-----------------------------------------------------------------------===*/

#ifndef TIMEVARYINGDRAGFORCE_H
#define TIMEVARYINGDRAGFORCE_H

#include "DragForce.h"

class TimeVaryingDragForce : public DragForce {
public:
	TimeVaryingDragForce(Cloud * const C, const double scale, const double offset)
	: DragForce(C, -offset), scaleConst(scale), offsetConst(offset) {}
	~TimeVaryingDragForce() {}

	void force1(const double currentTime); // rk substep 1
	void force2(const double currentTime); // rk substep 2
	void force3(const double currentTime); // rk substep 3
	void force4(const double currentTime); // rk substep 4
	
	void writeForce(fitsfile *file, int *error) const;
	void readForce(fitsfile *file, int *error);

private:
	double scaleConst, offsetConst; // [Hz/s], [Hz]

	const double calculateGamma(const double currentTime) const;
};

#endif // TIMEVARYINGDRAGFORCE_H
