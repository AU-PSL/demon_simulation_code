/*===- PositionVelocityCacheOperator.h - libSimulation -========================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef POSITIONVELOCITYCACHEOPERATOR_H
#define POSITIONVELOCITYCACHEOPERATOR_H

#include "Operator.h"

class PositionVelocityCacheOperator : public Operator {
public:
	PositionVelocityCacheOperator(Cloud * const mycloud) : Operator(mycloud) {}
	~PositionVelocityCacheOperator() {}
    
	void operation1(const double currentTime);
	void operation2(const double currentTime);
	void operation3(const double currentTime);
	void operation4(const double currentTime);
};

#endif // POSITIONVELOCITYCACHEOPERATOR
