/*===- ChargeOperator.h - libSimulation -=======================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef CHARGECACHEOPERATOR_H
#define CHARGECACHEOPERATOR_H

#include "Operator.h"

class ChargeOperator : public Operator 
{
public:
	ChargeOperator(Cloud * const mycloud) : Operator(mycloud) {}
	~ChargeOperator() {}
    
	void operation1(const double currentTime);
	void operation2(const double currentTime);
	void operation3(const double currentTime);
	void operation4(const double currentTime);
};

#endif // CHARGECACHEOPERATOR
