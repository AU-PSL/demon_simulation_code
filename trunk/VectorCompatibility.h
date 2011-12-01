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

#include <immintrin.h>

static inline void plusEqual_pd(double * const a, const __m128d b) {
	_mm_store_pd(a, _mm_load_pd(a) + b);
}

static inline void minusEqual_pd(double * const a, const __m128d b) {
	_mm_store_pd(a, _mm_load_pd(a) + b);
}

static inline void plusEqualr_pd(double * const a, const __m128d b) {
	_mm_storer_pd(a, _mm_loadr_pd(a) + b);
}

static inline void minusEqualr_pd(double * const a, const __m128d b) {
	_mm_storer_pd(a, _mm_loadr_pd(a) + b);
}

#if !defined(__GNUC__) && !defined(__clang__)
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
