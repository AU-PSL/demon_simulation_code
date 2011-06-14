/*===- FieldPotentialOperator.h - libSimulation -===============================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details. 
*
*===-----------------------------------------------------------------------===*/

#ifndef FIELDPOTENTIALOPERATOR_H
#define FIELDPOTENTIALOPERATOR_H

#include "Operator.h"
#include "VectorCompatibility.h"

class FieldPotentialOperator : public Operator
{	
public:
	FieldPotentialOperator(Cloud * const myCloud) : Operator(myCloud) {}
	~FieldPotentialOperator() {}

	void operation1(const double currentTime);
	void operation2(const double currentTime);
	void operation3(const double currentTime);
	void operation4(const double currentTime);

private:
	void potential(const cloud_index currentParticle, const double disp, const __m128d charges);
	void potential(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispV, const __m128d currentCharges, const __m128d iCharges);
	void potentialr(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispVr, const __m128d currentCharges, const __m128d iChargesr);

	void field(const cloud_index currentParticle, const __m128d dispV, const __m128d dispXv, const __m128d dispYv, const __m128d currentPhi);
	void field(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispV, const __m128d dispXv, const __m128d dispYv, const __m128d currentPhi, const __m128d iPhi);
	void fieldr(const cloud_index currentParticle, const cloud_index iParticle, const __m128d dispVr, const __m128d dispXvr, const __m128d dispYvr, const __m128d currentPhi, const __m128d iPhir);
};

#endif // FIELDPOTENTIALOPERATOR_H
