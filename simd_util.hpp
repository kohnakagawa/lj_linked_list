#pragma once

#include <x86intrin.h>
#include <iostream>

static inline void transpose_4x4(const __m256d& va,
                                 const __m256d& vb,
                                 const __m256d& vc,
                                 const __m256d& vd,
                                 __m256d& vx,
                                 __m256d& vy,
                                 __m256d& vz) {
  __m256d tmp0 = _mm256_unpacklo_pd(va, vb);
  __m256d tmp1 = _mm256_unpackhi_pd(va, vb);
  __m256d tmp2 = _mm256_unpacklo_pd(vc, vd);
  __m256d tmp3 = _mm256_unpackhi_pd(vc, vd);
  vx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
  vy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
  vz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
}

static void print256d(const __m256d& val) {
  union {
    __m256d v;
    double elem[4];
  } tmp;
  tmp.v = val;
  printf("%f %f %f %f\n",
         tmp.elem[0], tmp.elem[1], tmp.elem[2], tmp.elem[3]);
}
