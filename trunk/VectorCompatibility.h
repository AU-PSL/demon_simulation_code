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

#ifdef __AVX__
#define FLOAT_STRIDE  8
#define DOUBLE_STRIDE 4

typedef __m256d doubleV;
typedef __m256 floatV;

#else

#define FLOAT_STRIDE  4
#define DOUBLE_STRIDE 2

typedef __m128d doubleV;
typedef __m128 floatV;

#endif

/*===- Arithmatic ---------------------------------------------------------===*/

static inline void plusEqual_pd(double * const a, const doubleV b) {
#ifdef __AVX__
    _mm256_store_pd(a, _mm256_load_pd(a) + b);
#else
	_mm_store_pd(a, _mm_add_pd(_mm_load_pd(a), b));
#endif
}

static inline void minusEqual_pd(double * const a, const doubleV b) {
#ifdef __AVX__
    _mm256_store_pd(a, _mm256_load_pd(a) - b);
#else
	_mm_store_pd(a, _mm_sub_pd(_mm_load_pd(a), b));
#endif
}

static inline const doubleV fmadd_pd(const doubleV a, const doubleV b, const doubleV c) {
#ifdef __AVX__
    return _mm256_fmadd_pd(a, b, c);
#else
    return _mm_add_pd(_mm_mul_pd(a, b), c);
#endif
}

static inline const doubleV add_pd(const doubleV a, const doubleV b) {
#ifdef __AVX__
    return _mm256_add_pd(a, b);
#else
    return _mm_add_pd(a, b);
#endif
}

static inline const floatV add_ps(const floatV a, const floatV b) {
#ifdef __AVX__
    return _mm256_add_ps(a, b);
#else
    return _mm_add_ps(a, b);
#endif
}

static inline const doubleV sub_pd(const doubleV a, const doubleV b) {
#ifdef __AVX__
    return _mm256_sub_pd(a, b);
#else
    return _mm_sub_pd(a, b);
#endif
}

static inline const floatV sub_ps(const floatV a, const floatV b) {
#ifdef __AVX__
    return _mm256_sub_ps(a, b);
#else
    return _mm_sub_ps(a, b);
#endif
}

static inline const doubleV mul_pd(const doubleV a, const doubleV b) {
#ifdef __AVX__
    return _mm256_mul_pd(a, b);
#else
    return _mm_mul_pd(a, b);
#endif
}

static inline const floatV mul_ps(const floatV a, const floatV b) {
#ifdef __AVX__
    return _mm256_mul_ps(a, b);
#else
    return _mm_mul_ps(a, b);
#endif
}

static inline const doubleV div_pd(const doubleV a, const doubleV b) {
#ifdef __AVX__
    return _mm256_div_pd(a, b);
#else
    return _mm_div_pd(a, b);
#endif
}

static inline const doubleV sqrt_pd(const doubleV a) {
#ifdef __AVX__
    return _mm256_sqrt_pd(a);
#else
    return _mm_sqrt_pd(a);
#endif
}

static inline const floatV sqrt_ps(const floatV a) {
#ifdef __AVX__
    return _mm256_sqrt_ps(a);
#else
    return _mm_sqrt_ps(a);
#endif
}

/*===- Comparison ---------------------------------------------------------===*/

static inline const int movemask_ps(const floatV a) {
#ifdef __AVX__
    return _mm256_movemask_ps(a);
#else
    return _mm_movemask_ps(a);
#endif
}

static inline const floatV cmple_ps(const floatV a, const floatV b) {
#ifdef __AVX__
    return _mm256_cmp_ps(a, b, LE_OS);
#else
    return _mm_cmple_ps(a, b);
#endif
}

/*===- Comparison ---------------------------------------------------------===*/

static inline const doubleV set1_ps(const float a) {
#ifdef __AVX__
    return _mm256_set1_ps(a);
#else
    return _mm_set1_pd(a);
#endif
}

#endif // VECTORCOMPATIBILITY_H
