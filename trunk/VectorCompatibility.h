/*===- VectorCompatibility.h - libSimulation -==================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef VECTORCOMPATIBILITY_H
#define VECTORCOMPATIBILITY_H

#include <smmintrin.h>

#ifndef __GNUC__ || __clang__
__m128d operator+(const __m128d &a, const __m128d &b) {
	return _mm_add_pd(a, b);
}

__m128d operator-(const __m128d &a, const __m128d &b) {
	return _mm_sub_pd(a, b);
}

__m128d operator*(const __m128d &a, const __m128d &b) {
	return _mm_mul_pd(a, b);
}

__m128d operator/(const __m128d &a, const __m128d &b) {
	return _mm_div_pd(a, b);
}

__m128d operator&&(const __m128d &a, const __m128d &b) {
	return _mm_and_pd(a, b);
}
#endif

#endif // VECTORCOMPATIBILITY_H
