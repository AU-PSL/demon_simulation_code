/*===- DragForce.h - libSimulation -============================================
*
*                                  DEMON
* 
* This file is distributed under the BSD Open Source License. See LICENSE.TXT  
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef DRAGFORCE_H
#define DRAGFORCE_H

#include "Force.h"

class DragForce : public Force {
public:
	DragForce(Cloud * const myCloud, const double gamma)
	: Force(myCloud), dragConst(-gamma) {}
	~DragForce() {}

	virtual void force1(const double currentTime); // rk substep 1
	virtual void force2(const double currentTime); // rk substep 2
	virtual void force3(const double currentTime); // rk substep 3
	virtual void force4(const double currentTime); // rk substep 4
	
	virtual void writeForce(fitsfile * const file, int * const error) const;
	virtual void readForce(fitsfile * const file, int * const error);

protected:
	double dragConst; // [Hz]

private:
	void force(const cloud_index currentParticle, const doubleV currentVelocityX, const doubleV currentVelocityY);
};

#endif // DRAGFORCE_H
