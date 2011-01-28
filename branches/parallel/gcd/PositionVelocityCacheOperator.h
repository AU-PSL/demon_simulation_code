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
#include <dispatch/dispatch.h>

class PositionVelocityCacheOperator : public Operator 
{
public:
	PositionVelocityCacheOperator(Cloud * const mycloud) : Operator(mycloud), queue(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0)) {}
	~PositionVelocityCacheOperator() {}
    
	void operation1(const double currentTime);
	void operation2(const double currentTime);
	void operation3(const double currentTime);
	void operation4(const double currentTime);

private:
	dispatch_queue_t queue;
};

#endif // POSITIONVELOCITYCACHEOPERATOR
