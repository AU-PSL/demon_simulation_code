/*===- Operator.h - libSimulation -=============================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef OPERATOR_H
#define OPERATOR_H

#include "Cloud.h"

class Operator {
public:
	Cloud * const cloud;
	Operator(Cloud * const myCloud) : cloud(myCloud) {}
	virtual ~Operator() {}
    
	virtual void operation1(const double currentTime)=0;
	virtual void operation2(const double currentTime)=0;
	virtual void operation3(const double currentTime)=0;
	virtual void operation4(const double currentTime)=0;
};

#endif // OPERATOR_H
