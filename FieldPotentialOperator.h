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
	void fieldsAndPotentials(const cloud_index currentParticle, const __m128d currentCharge,
							 const double displacementX, const double displacementY);
	void fieldsAndPotentials(const cloud_index currentParticle, const cloud_index iParticle,
	                         const __m128d currentCharge, const __m128d iCharge,
	                         const __m128d displacementX, const __m128d displacementY);
	void fieldsAndPotentialsr(const cloud_index currentParticle, const cloud_index iParticle,
	                          const __m128d currentCharge, const __m128d iCharge,
	                          const __m128d displacementX, const __m128d displacementY);
	
	static const __m128d coulomb;
};

#endif // FIELDPOTENTIALOPERATOR_H
