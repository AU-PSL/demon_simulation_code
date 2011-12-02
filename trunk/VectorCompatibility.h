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

/*===- Load/Store ---------------------------------------------------------===*/

static inline void store_pd(double * const a, const doubleV b) {
#ifdef __AVX__
    return _mm256_store_pd(a, b);
#else
    return _mm_store_pd(a, b);
#endif
}

static inline const doubleV load_pd(double * const a) {
#ifdef __AVX__
    return _mm256_load_pd(a);
#else
    return _mm_load_pd(a);
#endif
}

/*===- Set ----------------------------------------------------------------===*/

static inline const doubleV set1_pd(const double a) {
#ifdef __AVX__
    return _mm256_set1_pd(a);
#else
    return _mm_set1_pd(a);
#endif
}

static inline const floatV set1_ps(const float a) {
#ifdef __AVX__
    return _mm256_set1_ps(a);
#else
    return _mm_set1_ps(a);
#endif
}

static inline const doubleV set0_pd() {
#ifdef __AVX__
    return _mm256_setzero_pd();
#else
    return _mm_setzero_pd();
#endif
}

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

static inline const doubleV sub_pd(const doubleV a, const double b) {
#ifdef __AVX__
    return _mm256_sub_pd(a, _mm256_set1_pd(b));
#else
    return _mm_sub_pd(a, _mm_set1_pd(b));
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

static inline const doubleV mul_pd(const doubleV a, const double b) {
#ifdef __AVX__
    return _mm256_mul_pd(a, _mm256_set1_pd(b));
#else
    return _mm_mul_pd(a, _mm_set1_pd(b));
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

static inline const doubleV div_pd(const doubleV a, const double b) {
#ifdef __AVX__
    return _mm256_div_pd(a, _mm256_set1_pd(b));
#else
    return _mm_div_pd(a, _mm_set1_pd(b));
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

static inline const int movemask_pd(const doubleV a) {
#ifdef __AVX__
    return _mm256_movemask_pd(a);
#else
    return _mm_movemask_pd(a);
#endif
}

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

static inline const doubleV cmpgt_pd(const doubleV a, const double b) {
#ifdef __AVX__
    return _mm256_cmp_pd(a, _mm256_set1_pd(b), LT_OS);
#else
    return _mm_cmpgt_pd(a, _mm_set1_pd(b));
#endif
}

static inline const doubleV cmplt_pd(const doubleV a, const double b) {
#ifdef __AVX__
    return _mm256_cmp_pd(a, _mm256_set1_pd(b), GT_OS);
#else
    return _mm_cmplt_pd(a, _mm_set1_pd(b));
#endif
}

/*===- Logical ------------------------------------------------------------===*/

static inline const doubleV and_pd(const doubleV a, const doubleV b) {
#ifdef __AVX__
    return _mm256_and_pd(a, b);
#else
    return _mm_and_pd(a, b);
#endif
}

/*===- math functions -----------------------------------------------------===*/

static inline const doubleV exp_pd(const doubleV a) {
    double b[DOUBLE_STRIDE];
    store_pd(b, a);
    
#ifdef __AVX__
    return _mm256_set_pd(exp(b[3]), exp(b[2]), exp(b[1]), exp(b[0]));
#else
    return _mm_set_pd(exp(b[1]), exp(b[0]));
#endif
}

static inline const doubleV sin_pd(const doubleV a) {
    double b[DOUBLE_STRIDE];
    store_pd(b, a);
    
#ifdef __AVX__
    return _mm256_set_pd(sin(b[3]), sin(b[2]), sin(b[1]), sin(b[0]));
#else
    return _mm_set_pd(sin(b[1]), sin(b[0]));
#endif
}

static inline const doubleV cos_pd(const doubleV a) {
    double b[DOUBLE_STRIDE];
    store_pd(b, a);
    
#ifdef __AVX__
    return _mm256_set_pd(cos(b[3]), cos(b[2]), cos(b[1]), cos(b[0]));
#else
    return _mm_set_pd(cos(b[1]), cos(b[0]));
#endif
}

static inline const doubleV length_pd(const doubleV a, const doubleV b) {
    return sqrt_pd(add_pd(mul_pd(a, a), mul_pd(b, b)));
}

/*===- Misc ---------------------------------------------------------------===*/

static inline const doubleV select_pd(const int mask, const double trueValue, const double falseValue) {
#ifdef __AVX__
    return _mm256_set_pd((mask & 8) ? trueValue : falseValue, 
                         (mask & 4) ? trueValue : falseValue,
                         (mask & 2) ? trueValue : falseValue, // _mm256_set_pd() is backwards.
                         (mask & 1) ? trueValue : falseValue);
#else
    return _mm_set_pd((mask & 2) ? trueValue : falseValue, // _mm_set_pd() is backwards.
                      (mask & 1) ? trueValue : falseValue);
#endif
}

#endif // VECTORCOMPATIBILITY_H
