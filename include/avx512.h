
/*******************************************************************************
*
* File avx512.h
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Macros for operating on SU(3) vectors and matrices using Intel intrinsic
* operations for AVX512 compatible processors
*
*******************************************************************************/

#ifndef AVX512_H
#define AVX512_H

#ifndef SSE2_H
#include "sse2.h"
#endif

#include "immintrin.h"



/* Macros for single precision floating point numbers */

/* Write 6 color weyl vectors as a spinor */
#define _avx512_write_6_hwv_f( c1,c2,c3, a ){            \
  __m256 t256; __m128 t128;                              \
  t256 = _mm256_shuffle_ps( c1, c2, 0b01000100 );        \
  t128 = _mm256_castps256_ps128( t256 );                 \
  _mm_storeu_ps( a, t128 );                               \
  t128 = _mm256_extractf128_ps( t256, 1 );               \
  _mm_storeu_ps( a+12, t128 );                            \
                                                         \
  t256 = _mm256_shuffle_ps( c3, c1, 0b11100100 );        \
  t128 = _mm256_castps256_ps128( t256 );                 \
  _mm_storeu_ps( a+4, t128 );                             \
  t128 = _mm256_extractf128_ps( t256, 1 );               \
  _mm_storeu_ps( a+16, t128 );                            \
                                                         \
  t256 = _mm256_shuffle_ps( c2, c3, 0b11101110 );        \
  t128 = _mm256_castps256_ps128( t256 );                 \
  _mm_storeu_ps( a+8, t128 );                             \
  t128 = _mm256_extractf128_ps( t256, 1 );               \
  _mm_storeu_ps( a+20, t128 );                            \
}


/* Load 6 color weyl vectors from a spinor */
#define _avx512_load_6_hwv_f( c1,c2,c3, a ){             \
  __m256 t1,t2,t3; __m128 t11, t14;                      \
  t11 = _mm_loadu_ps( a );                               \
  t14 = _mm_loadu_ps( a+12 );                            \
  t1 = _mm256_castps128_ps256( t11 );                    \
  t1 = _mm256_insertf128_ps( t1, t14, 1 );               \
                                                         \
  t11 = _mm_loadu_ps( a+4 );                             \
  t14 = _mm_loadu_ps( a+16 );                            \
  t2 = _mm256_castps128_ps256( t11 );                    \
  t2 = _mm256_insertf128_ps( t2, t14, 1 );               \
                                                         \
  t11 = _mm_loadu_ps( a+8 );                             \
  t14 = _mm_loadu_ps( a+20 );                            \
  t3 = _mm256_castps128_ps256( t11 );                    \
  t3 = _mm256_insertf128_ps( t3, t14, 1 );               \
                                                         \
  c1 = _mm256_shuffle_ps( t1, t2, 0b11100100 );          \
  c2 = _mm256_shuffle_ps( t1, t3, 0b01001110 );          \
  c3 = _mm256_shuffle_ps( t2, t3, 0b11100100 );          \
}



/* Load 4x2 complex numbers into an aray */
#define _avx512_load_4_cf( r, c1,c2,c3,c4 ){            \
     __m128 t128;                                       \
     t128 = _mm_loadu_ps( c1 );                         \
     r = _mm512_castps128_ps512( t128 );                \
     t128 = _mm_loadu_ps( c2 );                         \
     r = _mm512_insertf32x4 ( r, t128, 1);              \
     t128 = _mm_loadu_ps( c3 );                         \
     r = _mm512_insertf32x4 ( r, t128, 2);              \
     t128 = _mm_loadu_ps( c4 );                         \
     r = _mm512_insertf32x4 ( r, t128, 3);              \
}

/* Load 4 half-spinors and organize colorwise into vectors r1, r2 and r3 */
#define _avx512_load_4_halfspinor_f( r1, r2, r3, s1, s2, s3, s4   )   \
{                                                                     \
     __m512 t512a, t512b, t512c;                                      \
     _avx512_load_4_cf( t512a, s1,s2,s3,s4 );                         \
     _avx512_load_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                 \
     _avx512_load_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                 \
                                                                      \
     r1 = _mm512_shuffle_ps(t512a,t512b, 0b11100100);                 \
     r2 = _mm512_shuffle_ps(t512a,t512c, 0b01001110);                 \
     r3 = _mm512_shuffle_ps(t512b,t512c, 0b11100100);                 \
}

/* Load 4 half-spinors reversing the second two spinor indeces and
 * organize colorwise into vectors r1, r2 and r3 */
#define _avx512_load_4_halfspinor_f_reverse_up( r1, r2, r3, s1, s2, s3, s4 )  \
{                                                                             \
     __m512 t512a, t512b, t512c;                                              \
     __m512i idx;                                                             \
     _avx512_load_4_cf( t512a, s1,s2,s3,s4 );                                 \
     _avx512_load_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                         \
     _avx512_load_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                         \
                                                                              \
     idx = _mm512_setr_epi32( 0,1,16+2,16+3, 4,5,16+6,16+7,                   \
                              16+10,16+11,8,9, 16+14,16+15,12,13  );          \
     r1 = _mm512_permutex2var_ps( t512a, idx, t512b );                        \
     idx = _mm512_setr_epi32( 2,3,16+0,16+1, 6,7,16+4,16+5,                   \
                              16+8,16+9,10,11, 16+12,16+13,14,15  );          \
     r2 = _mm512_permutex2var_ps( t512a, idx, t512c );                        \
     idx = _mm512_setr_epi32( 0,1,16+2,16+3, 4,5,16+6,16+7,                   \
                              16+10,16+11,8,9, 16+14,16+15,12,13  );          \
     r3 = _mm512_permutex2var_ps( t512b, idx, t512c );                        \
}

/* Load 4 half-spinors reversing first two the spinor indeces and
 * organize colorwise into vectors r1, r2 and r3 */
#define _avx512_load_4_halfspinor_f_reverse_dn( r1, r2, r3, s1, s2, s3, s4 )  \
{                                                                             \
     __m512 t512a, t512b, t512c;                                              \
     __m512i idx;                                                             \
     _avx512_load_4_cf( t512a, s1,s2,s3,s4 );                                 \
     _avx512_load_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                         \
     _avx512_load_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                         \
                                                                              \
     idx = _mm512_setr_epi32( 16+2,16+3,0,1, 16+6,16+7,4,5,                   \
                              8,9,16+10,16+11, 12,13,16+14,16+15  );          \
     r1 = _mm512_permutex2var_ps( t512a, idx, t512b );                        \
     idx = _mm512_setr_epi32( 16+0,16+1,2,3, 16+4,16+5,6,7,                   \
                              10,11,16+8,16+9, 14,15,16+12,16+13  );          \
     r2 = _mm512_permutex2var_ps( t512a, idx, t512c );                        \
     idx = _mm512_setr_epi32( 16+2,16+3,0,1, 16+6,16+7,4,5,                   \
                              8,9,16+10,16+11, 12,13,16+14,16+15  );          \
     r3 = _mm512_permutex2var_ps( t512b, idx, t512c );                        \
}

/* Load 4x2 complex numbers into an aray */
#define _avx512_write_4_cf( r, c1,c2,c3,c4 ){                          \
     __m512 t512;                                                      \
     __m128 t128 = _mm512_castps512_ps128( r );                        \
     _mm_storeu_ps( c1, t128 );                                        \
     t128 = _mm512_extractf32x4_ps( r, 1 );                            \
     _mm_storeu_ps( c2, t128 );                                        \
     t128 = _mm512_extractf32x4_ps( r, 2 );                            \
     _mm_storeu_ps( c3, t128 );                                        \
     t128 = _mm512_extractf32x4_ps( r, 3 );                            \
     _mm_storeu_ps( c4, t128 );                                        \
}

/* Store 4 half-spinors from color vectors */
#define _avx512_write_4_halfspinor_f( r1, r2, r3, s1, s2, s3, s4   )   \
{                                                                      \
     __m512 t512a, t512b, t512c;                                       \
                                                                       \
     t512a = _mm512_shuffle_ps(r1,r2, 0b01000100);                     \
     t512b = _mm512_shuffle_ps(r3,r1, 0b11100100);                     \
     t512c = _mm512_shuffle_ps(r2,r3, 0b11101110);                     \
                                                                       \
     _avx512_write_4_cf( t512a, s1,s2,s3,s4 );                         \
     _avx512_write_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                 \
     _avx512_write_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                 \
}

/* Store 4 half-spinors from color vectors reversing the first two Dirac indeces */
#define _avx512_write_4_halfspinor_f_reverse_up( r1, r2, r3, s1, s2, s3, s4 ) \
{                                                                             \
     __m512 t512a, t512b, t512c;                                              \
     __m512i idx;                                                             \
     idx = _mm512_setr_epi32( 0,1,16+0,16+1, 4,5,16+4,16+5,                   \
                              10,11,16+10,16+11, 14,15,16+14,16+15  );        \
     t512a = _mm512_permutex2var_ps( r1, idx, r2 );                           \
     idx = _mm512_setr_epi32( 0,1,16+2,16+3, 4,5,16+6,16+7,                   \
                              10,11,16+8,16+9, 14,15,16+12,16+13  );          \
     t512b = _mm512_permutex2var_ps( r3, idx, r1 );                           \
     idx = _mm512_setr_epi32( 2,3,16+2,16+3, 6,7,16+6,16+7,                   \
                              8,9,16+8,16+9, 12,13,16+12,16+13  );            \
     t512c = _mm512_permutex2var_ps( r2, idx, r3 );                           \
                                                                              \
     _avx512_write_4_cf( t512a, s1,s2,s3,s4 );                                \
     _avx512_write_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                        \
     _avx512_write_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                        \
}


/* Store 4 half-spinors from color vectors reversing the second two Dirac indeces */
#define _avx512_write_4_halfspinor_f_reverse_dn( r1, r2, r3, s1, s2, s3, s4 ) \
{                                                                             \
     __m512 t512a, t512b, t512c;                                              \
     __m512i idx;                                                             \
     idx = _mm512_setr_epi32( 2,3,16+2,16+3, 6,7,16+6,16+7,                   \
                              8,9,16+8,16+9, 12,13,16+12,16+13  );            \
     t512a = _mm512_permutex2var_ps( r1, idx, r2 );                           \
     idx = _mm512_setr_epi32( 2,3,16+0,16+1, 6,7,16+4,16+5,                   \
                              8,9,16+10,16+11, 12,13,16+14,16+15  );          \
     t512b = _mm512_permutex2var_ps( r3, idx, r1 );                           \
     idx = _mm512_setr_epi32( 0,1,16+0,16+1, 4,5,16+4,16+5,                   \
                              10,11,16+10,16+11, 14,15,16+14,16+15  );        \
     t512c = _mm512_permutex2var_ps( r2, idx, r3 );                           \
                                                                              \
     _avx512_write_4_cf( t512a, s1,s2,s3,s4 );                                \
     _avx512_write_4_cf( t512b, s1+4,s2+4,s3+4,s4+4 );                        \
     _avx512_write_4_cf( t512c, s1+8,s2+8,s3+8,s4+8 );                        \
}








/* Combine Dirac spinors to half-spinors in deo and doe.
 */
#define _avx512_dirac_combine_f_1( a, b )                                \
{                                                                        \
  __m512i indexes;                                                       \
  __m512 c;                                                              \
  indexes = _mm512_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7,                   \
                               9, 8, 11, 10, 13, 12, 15, 14 );           \
  c = _mm512_permutexvar_ps( indexes, b);                                \
  a = _mm512_mask_add_ps( a, 0b0101101000001111, a, c );                 \
  a = _mm512_mask_sub_ps( a, 0b1010010111110000, a, c );                 \
}

#define _avx512_dirac_combine_f_2( a, b )                                \
{                                                                        \
  __m512i indexes;                                                       \
  __m512 c;                                                              \
  indexes = _mm512_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7,                   \
                               9, 8, 11, 10, 13, 12, 15, 14 );           \
  c = _mm512_permutexvar_ps( indexes, b);                                \
  a = _mm512_mask_add_ps( a, 0b1001011011000011, a, c );                 \
  a = _mm512_mask_sub_ps( a, 0b0110100100111100, a, c );                 \
}


#define _avx512_dirac_combine_f_3( a, b )                                \
{                                                                        \
  __m512i indexes;                                                       \
  __m512 c;                                                              \
  indexes = _mm512_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7,                   \
                               9, 8, 11, 10, 13, 12, 15, 14 );           \
  c = _mm512_permutexvar_ps( indexes, b);                                \
  a = _mm512_mask_add_ps( a, 0b1010010100001111, a, c );                 \
  a = _mm512_mask_sub_ps( a, 0b0101101011110000, a, c );                 \
}


#define _avx512_dirac_combine_f_4( a, b )                                \
{                                                                        \
  __m512i indexes;                                                       \
  __m512 c;                                                              \
  indexes = _mm512_setr_epi32( 0, 1, 2, 3, 4, 5, 6, 7,                   \
                       9, 8, 11, 10, 13, 12, 15, 14 );                   \
  c = _mm512_permutexvar_ps( indexes, b);                                \
  a = _mm512_mask_add_ps( a, 0b0110100111000011, a, c );                 \
  a = _mm512_mask_sub_ps( a, 0b1001011000111100, a, c );                 \
}




/* Multiply 4 vectors with a su(3) matrices, taking the inverse of every
 * second matrix
 */
#define avx512_su3_mixed_multiply_8( u1, um1, u2, um2, b1,b2,b3, a1,a2,a3 )           \
 {                                                                                    \
   __m512 ut11, ut21, ut31, ut41, ut12, ut22, ut32, ut42;                             \
   __m512 ut1, ut2, sign, c;                                                          \
   __m512 t1,t2,t3;                                                                   \
   __m512i indexes;                                                                   \
   ut11 = _mm512_loadu_ps( &(u1).c11.re );  ut12 = _mm512_loadu_ps( &(u1).c33.re );   \
   ut21 = _mm512_loadu_ps( &(um1).c11.re ); ut22 = _mm512_loadu_ps( &(um1).c33.re );  \
   ut31 = _mm512_loadu_ps( &(u2).c11.re );  ut32 = _mm512_loadu_ps( &(u2).c33.re );   \
   ut41 = _mm512_loadu_ps( &(um2).c11.re ); ut42 = _mm512_loadu_ps( &(um2).c33.re );  \
                                                                                      \
   indexes = _mm512_setr_epi32( 0, 1, 2, 3, 8, 9, 6, 7,                               \
                                16, 17, 18, 19, 24, 25, 22, 23 );                     \
   ut1 = _mm512_permutex2var_ps( ut11, indexes, ut21 );                               \
   ut2 = _mm512_permutex2var_ps( ut31, indexes, ut41 );                               \
                                                                                      \
   indexes = _mm512_setr_epi32( 0,0,0,0, 8,8,8,8, 16,16,16,16, 24,24,24,24 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_mul_ps ( a1, c );                                                      \
                                                                                      \
   indexes = _mm512_setr_epi32( 4,4,4,4, 12,12,12,12, 20,20,20,20, 28,28,28,28 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_mul_ps ( a2, c );                                                      \
                                                                                      \
   indexes = _mm512_setr_epi32( 2,2,2,2, 14,14,14,14, 18,18,18,18, 30,30,30,30 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_fmadd_ps ( a2, c, b1 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 6,6,6,6, 10,10,10,10, 22,22,22,22, 26,26,26,26 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_fmadd_ps ( a1, c, b2 );                                                \
                                                                                      \
                                                                                      \
   sign = _mm512_set_ps( -1,1,-1,1, 1,-1,1,-1, -1,1,-1,1, 1,-1,1,-1 );                \
   t1 = _mm512_permute_ps( a1, 0b10110001  );                                         \
   t1 = _mm512_mul_ps( t1, sign );                                                    \
   t2 = _mm512_permute_ps( a2, 0b10110001  );                                         \
   t2 = _mm512_mul_ps( t2, sign );                                                    \
                                                                                      \
   indexes = _mm512_setr_epi32( 1,1,1,1, 9,9,9,9, 17,17,17,17, 25,25,25,25 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_fmadd_ps ( t1, c, b1 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 5,5,5,5, 13,13,13,13, 21,21,21,21, 29,29,29,29 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_fmadd_ps ( t2, c, b2 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 3,3,3,3, 15,15,15,15, 19,19,19,19, 31,31,31,31 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_fmadd_ps ( t2, c, b1 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 7,7,7,7, 11,11,11,11, 23,23,23,23, 27,27,27,27 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_fmadd_ps ( t1, c, b2 );                                                \
                                                                                      \
                                                                                      \
   indexes = _mm512_setr_epi32( 4, 5, 12, 13, 10, 11, 14, 15,                         \
                                20, 21, 28, 29, 26, 27, 30, 31 );                     \
   ut1 = _mm512_permutex2var_ps( ut11, indexes, ut21 );                               \
   ut2 = _mm512_permutex2var_ps( ut31, indexes, ut41 );                               \
                                                                                      \
   indexes = _mm512_setr_epi32( 0,0,0,0, 10,10,10,10, 16,16,16,16, 26,26,26,26 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_fmadd_ps ( a3, c, b1 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 2,2,2,2, 8,8,8,8, 18,18,18,18, 24,24,24,24 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_mul_ps ( a1, c );                                                      \
                                                                                      \
   indexes = _mm512_setr_epi32( 4,4,4,4, 14,14,14,14, 20,20,20,20, 30,30,30,30 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_fmadd_ps ( a3, c, b2 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 6,6,6,6, 12,12,12,12, 22,22,22,22, 28,28,28,28 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_fmadd_ps ( a2, c, b3 );                                                \
                                                                                      \
                                                                                      \
   t3 = _mm512_permute_ps( a3, 0b10110001  );                                         \
   t3 = _mm512_mul_ps( t3, sign );                                                    \
                                                                                      \
   indexes = _mm512_setr_epi32( 1,1,1,1, 11,11,11,11, 17,17,17,17, 27,27,27,27 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b1 = _mm512_fmadd_ps ( t3, c, b1 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 3,3,3,3, 9,9,9,9, 19,19,19,19, 25,25,25,25 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_fmadd_ps ( t1, c, b3 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 5,5,5,5, 15,15,15,15, 21,21,21,21, 31,31,31,31 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b2 = _mm512_fmadd_ps ( t3, c, b2 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 7,7,7,7, 13,13,13,13, 23,23,23,23, 29,29,29,29 );     \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_fmadd_ps ( t2, c, b3 );                                                \
                                                                                      \
                                                                                      \
   indexes = _mm512_setr_epi32( 0, 1, 16, 17, 0,0,0,0, 0,0,0,0, 0,0,0,0 );            \
   ut1 = _mm512_permutex2var_ps( ut12, indexes, ut22 );                               \
   ut2 = _mm512_permutex2var_ps( ut32, indexes, ut42 );                               \
                                                                                      \
   indexes = _mm512_setr_epi32( 0,0,0,0, 2,2,2,2, 16,16,16,16, 18,18,18,18 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_fmadd_ps ( a3, c, b3 );                                                \
                                                                                      \
   indexes = _mm512_setr_epi32( 1,1,1,1, 3,3,3,3, 17,17,17,17, 19,19,19,19 );         \
   c = _mm512_permutex2var_ps( ut1, indexes, ut2 );                                   \
   b3 = _mm512_fmadd_ps ( t3, c, b3 );                                                \
}



/* Insert 256 bits into a 512 bit single precision vector */
#define _avx512_insert_256_h_f( a, t )                             \
{                                                                  \
  __m512d td512;                                                   \
  __m256d td256;                                                   \
  td512 = _mm512_castps_pd( a );                                   \
  td256 = _mm256_castps_pd( t );                                   \
  td512 = _mm512_insertf64x4( td512, td256, 1 );                   \
  a = _mm512_castpd_ps( td512 );                                   \
}

/* Extract 256 bits from a 512 bit single precision vector */
#define _avx512_extract_256_h_f( t, a )                           \
{                                                                 \
  __m512d td512;                                                  \
  __m256d td256;                                                  \
  td512 = _mm512_castps_pd( a );                                  \
  td256 = _mm256_castps_pd( t );                                  \
  td256 = _mm512_extractf64x4_pd( td512, 1 );                     \
  t = _mm256_castpd_ps( td256 );                                  \
}



/* Accumulate elements of Dirac vectors into a Weyl vector in deo and doe */
#define _avx512_to_weyl_f_12( c, b ){                                  \
  __m256 w, sign, t5,t6;                                               \
  __m256i idx;                                                         \
  t5 = _mm512_castps512_ps256( b );                                    \
                                                                       \
  idx = _mm256_setr_epi32( 4, 5, 6, 7, 0, 1, 2, 3 );                   \
  t6 = _mm256_permutevar8x32_ps( t5, idx );                            \
  sign = _mm256_set_ps( -1,-1,-1,-1, 1,1,1,1 );                        \
  c = _mm256_fmadd_ps( t5, sign, t6 );                                 \
                                                                       \
  _avx512_extract_256_h_f( t5, b );                                    \
  t6 = _mm256_permutevar8x32_ps( t5, idx );                            \
  sign = _mm256_set_ps( -1,-1,-1,-1, 1,1,1,1 );                        \
  w = _mm256_fmadd_ps( t6, sign, t5 );                                 \
  idx = _mm256_setr_epi32( 0, 1, 2, 3,  3, 2, 1, 0 );                  \
  w = _mm256_permutevar_ps( w, idx );                                  \
  sign = _mm256_set_ps( 1,-1,1,-1, 1,1,1,1 );                          \
  c = _mm256_fmadd_ps( w, sign, c );                                   \
}

#define _avx512_to_weyl_f_34( c, b ){                                  \
  __m256 w, sign, t5,t6;                                               \
  __m256i idx;                                                         \
  t5 = _mm512_castps512_ps256( b );                                    \
                                                                       \
  idx = _mm256_setr_epi32( 4, 5, 6, 7, 0, 1, 2, 3 );                   \
  t6 = _mm256_permutevar8x32_ps( t5, idx );                            \
  sign = _mm256_set_ps( -1,-1,-1,-1, 1,1,1,1 );                        \
  w = _mm256_fmadd_ps( t6, sign, t5 );                                 \
  idx = _mm256_setr_epi32( 0, 1, 2, 3,  2, 3, 0, 1 );                  \
  w = _mm256_permutevar_ps( w, idx );                                  \
  sign = _mm256_set_ps( -1,-1,1,1, 1,1,1,1 );                          \
  c = _mm256_fmadd_ps( w, sign, c );                                   \
                                                                       \
  _avx512_extract_256_h_f( t5, b );                                    \
  idx = _mm256_setr_epi32( 4, 5, 6, 7, 0, 1, 2, 3 );                   \
  t6 = _mm256_permutevar8x32_ps( t5, idx );                            \
  sign = _mm256_set_ps( -1,-1,-1,-1, 1,1,1,1 );                        \
  w = _mm256_fmadd_ps( t6, sign, t5 );                                 \
  idx = _mm256_setr_epi32( 0, 1, 2, 3,  1, 0, 3, 2 );                  \
  w = _mm256_permutevar_ps( w, idx );                                  \
  sign = _mm256_set_ps( -1,1,1,-1, 1,1,1,1 );                          \
  c = _mm256_fmadd_ps( w, sign, c );                                   \
}


/* Expand a Weyl vector into a Dirac vector in deo and doe */
#define _avx512_to_dirac_f_1( a1, c1 )                                 \
{                                                                      \
  __m256 t1,t2, sign;                                                  \
  __m256i idx;                                                         \
  idx = _mm256_setr_epi32( 4, 5, 6, 7, 0, 1, 2, 3 );                   \
  t1 = _mm256_permutevar8x32_ps( c1, idx );                            \
  sign = _mm256_set_ps( -1,-1,-1,-1, 1,1,1,1 );                        \
  t2 = _mm256_fmadd_ps( c1, sign, t1 );                                \
  a1 = _mm512_castps256_ps512( t2 );                                   \
}


#define _avx512_to_dirac_f_2( a1, c1 )                                 \
{                                                                      \
  __m256 t1,t2, sign;                                                  \
  __m256i idx;                                                         \
  idx = _mm256_setr_epi32( 7, 6, 5, 4, 7, 6, 5, 4 );                   \
  t1 = _mm256_permutevar8x32_ps( c1, idx );                            \
  idx = _mm256_setr_epi32( 0, 1, 2, 3, 0, 1, 2, 3 );                   \
  t2 = _mm256_permutevar8x32_ps( c1, idx );                            \
  sign = _mm256_set_ps( -1,1,-1,1, 1,-1,1,-1 );                        \
  t2 = _mm256_fmadd_ps( t1, sign, t2 );                                \
  _avx512_insert_256_h_f( a1, t2 );                                    \
}

#define _avx512_to_dirac_f_3( a1, c1 )                                 \
{                                                                      \
  __m512 t5,t6,t7;                                                     \
  __m512i idx;                                                         \
  t5 = _mm512_castps256_ps512( c1 );                                   \
  idx = _mm512_setr_epi32( 6,7,4,5, 6,7,4,5, 5,4,7,6, 5,4,7,6 );       \
  t6 = _mm512_permutexvar_ps( idx, t5 );                               \
  idx = _mm512_setr_epi32( 0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3 );       \
  t7 = _mm512_permutexvar_ps( idx, t5 );                               \
  a1 = _mm512_maskz_add_ps(    0b1001011011000011, t7, t6 );           \
  a1 = _mm512_mask_sub_ps( a1, 0b0110100100111100, t7, t6 );           \
}















/* Macros for double precision numbers */

/* Load 2 half-spinors and organize colorwise into vectors s1, s2 and s3 */
#define _avx512_load_2_halfspinor_d( s1, s2, s3, sp, sm )            \
{                                                                    \
  __m512d a1,a2,a3,a4,a5,a6;                                         \
  __m512i idx;                                                       \
  a1 = _mm512_loadu_pd( sm );                                        \
  a2 = _mm512_loadu_pd( sm+8 );                                      \
  a3 = _mm512_loadu_pd( sp );                                        \
  a4 = _mm512_loadu_pd( sp+8 );                                      \
                                                                     \
  idx = _mm512_setr_epi64( 8,9,14,15, 0,1,6,7 );                     \
  s1 = _mm512_permutex2var_pd( a1, idx, a3 );                        \
                                                                     \
  idx = _mm512_setr_epi64( 4,5,10,11, 2,3,8,9 );                     \
  a5 = _mm512_permutex2var_pd( a1, idx, a2 );                        \
  idx = _mm512_setr_epi64( 2,3,8,9, 4,5,10,11 );                     \
  a6 = _mm512_permutex2var_pd( a3, idx, a4 );                        \
  idx = _mm512_setr_epi64( 0,1,2,3, 12,13,14,15 );                   \
  s2 = _mm512_permutex2var_pd( a6, idx, a5 );                        \
                                                                     \
  idx = _mm512_setr_epi64( 4,5,6,7, 8,9,10,11 );                     \
  s3 = _mm512_permutex2var_pd( a6, idx, a5 );                        \
}

/* Load 2 half-spinors reversing the spinor indeces and
 * organize colorwise into vectors s1, s2 and s3 */
#define _avx512_load_2_halfspinor_d_reverse( s1, s2, s3, sp, sm )    \
{                                                                    \
  __m512d a1,a2,a3,a4,a5,a6;                                         \
  __m512i idx;                                                       \
  a1 = _mm512_loadu_pd( sm );                                        \
  a2 = _mm512_loadu_pd( sm+8 );                                      \
  a3 = _mm512_loadu_pd( sp );                                        \
  a4 = _mm512_loadu_pd( sp+8 );                                      \
                                                                     \
  idx = _mm512_setr_epi64( 14,15,8,9, 6,7,0,1 );                     \
  s1 = _mm512_permutex2var_pd( a1, idx, a3 );                        \
                                                                     \
  idx = _mm512_setr_epi64( 10,11,4,5, 8,9,2,3 );                     \
  a5 = _mm512_permutex2var_pd( a1, idx, a2 );                        \
  idx = _mm512_setr_epi64( 8,9,2,3, 10,11,4,5 );                     \
  a6 = _mm512_permutex2var_pd( a3, idx, a4 );                        \
  idx = _mm512_setr_epi64( 0,1,2,3, 12,13,14,15 );                   \
  s2 = _mm512_permutex2var_pd( a6, idx, a5 );                        \
                                                                     \
  idx = _mm512_setr_epi64( 4,5,6,7, 8,9,10,11 );                     \
  s3 = _mm512_permutex2var_pd( a6, idx, a5 );                        \
}

/* Write 2 half-spinors from three color vectors */
#define _avx512_store_2_halfspinor_d( s1, s2, s3, sp, sm )           \
{                                                                    \
  __m512d a1,a2,a3,a4,a5,a6;                                         \
  __m256d l;                                                         \
  __m512i idx;                                                       \
  idx = _mm512_setr_epi64( 0,1,8,9, 4,5,12,13 );                     \
  a1 = _mm512_permutex2var_pd( s1, idx, s2 );                        \
  idx = _mm512_setr_epi64( 0,1,10,11, 4,5,14,15 );                   \
  a2 = _mm512_permutex2var_pd( s3, idx, s1 );                        \
  idx = _mm512_setr_epi64( 2,3,10,11, 6,7,14,15 );                   \
  a3 = _mm512_permutex2var_pd( s2, idx, s3 );                        \
                                                                     \
  l = _mm512_castpd512_pd256( a1 );                                  \
  _mm256_storeu_pd( sp, l );                                          \
  l = _mm512_castpd512_pd256( a2 );                                  \
  _mm256_storeu_pd( sp+4, l );                                        \
  l = _mm512_castpd512_pd256( a3 );                                  \
  _mm256_storeu_pd( sp+8, l );                                        \
                                                                     \
  l = _mm512_extractf64x4_pd( a1, 1 );                               \
  _mm256_storeu_pd( sm, l );                                          \
  l = _mm512_extractf64x4_pd( a2, 1 );                               \
  _mm256_storeu_pd( sm+4, l );                                        \
  l = _mm512_extractf64x4_pd( a3, 1 );                               \
  _mm256_storeu_pd( sm+8, l );                                        \
}



/* Multiply the lower half of the color vectors distributed in c1, c2 and c3
 * by the su3 matrix u and the upper half by the conjugate of um
 * Store in b1, b2 and b3
 */
#define avx512_su3_mul_quad_dble( u, um, b1, b2, b3, c1, c2, c3   )     \
{                                                                       \
  __m512d tu1, tu2, tu3, tum1, tum2, tum3;                              \
  __m512d u1;                              \
  __m512d t1, t2, t3, sign;                                             \
  __m512i indexes;                                                      \
  tu1 = _mm512_loadu_pd( &(u).c11.re );                                 \
  tu2 = _mm512_loadu_pd( &(u).c22.re );                                 \
  tu3 = _mm512_loadu_pd( &(u).c33.re );                                 \
  tum1 = _mm512_loadu_pd( &(um).c11.re );                               \
  tum2 = _mm512_loadu_pd( &(um).c22.re );                               \
  tum3 = _mm512_loadu_pd( &(um).c33.re );                               \
                                                                        \
  sign = _mm512_set_pd( -1,1,-1,1, 1,-1,1,-1 );                         \
  t1 = _mm512_permute_pd( c1, 0b01010101 );                             \
  t2 = _mm512_permute_pd( c2, 0b01010101 );                             \
  t3 = _mm512_permute_pd( c3, 0b01010101 );                             \
  t1 = _mm512_mul_pd( t1, sign );                                       \
  t2 = _mm512_mul_pd( t2, sign );                                       \
  t3 = _mm512_mul_pd( t3, sign );                                       \
                                                                        \
  indexes = _mm512_setr_epi64( 0, 0, 0, 0, 8, 8, 8, 8 );                \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b1 = _mm512_mul_pd ( u1, c1 );                                        \
  indexes = _mm512_setr_epi64( 1, 1, 1, 1, 9, 9, 9, 9 );                \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b1 = _mm512_fmadd_pd ( u1, t1, b1 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 2, 2, 2, 2, 14, 14, 14, 14 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b1 = _mm512_fmadd_pd ( u1, c2, b1 );                                  \
  indexes = _mm512_setr_epi64( 3, 3, 3, 3, 15, 15, 15, 15 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b1 = _mm512_fmadd_pd ( u1, t2, b1 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 4, 4, 4, 4, 12, 12, 12, 12 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum2 );                    \
  b1 = _mm512_fmadd_pd ( u1, c3, b1 );                                  \
  indexes = _mm512_setr_epi64( 5, 5, 5, 5, 13, 13, 13, 13 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum2 );                    \
  b1 = _mm512_fmadd_pd ( u1, t3, b1 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 6, 6, 6, 6, 10, 10, 10, 10 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b2 = _mm512_mul_pd ( u1, c1 );                                        \
  indexes = _mm512_setr_epi64( 7, 7, 7, 7, 11, 11, 11, 11 );            \
  u1 = _mm512_permutex2var_pd( tu1, indexes, tum1 );                    \
  b2 = _mm512_fmadd_pd ( u1, t1, b2 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 0, 0, 0, 0, 8, 8, 8, 8 );                \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b2 = _mm512_fmadd_pd ( u1, c2, b2 );                                  \
  indexes = _mm512_setr_epi64( 1, 1, 1, 1, 9, 9, 9, 9 );                \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b2 = _mm512_fmadd_pd ( u1, t2, b2 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 2, 2, 2, 2, 14, 14, 14, 14 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b2 = _mm512_fmadd_pd ( u1, c3, b2 );                                  \
  indexes = _mm512_setr_epi64( 3, 3, 3, 3, 15, 15, 15, 15 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b2 = _mm512_fmadd_pd ( u1, t3, b2 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 4, 4, 4, 4, 12, 12, 12, 12 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum1 );                    \
  b3 = _mm512_mul_pd ( u1, c1 );                                        \
  indexes = _mm512_setr_epi64( 5, 5, 5, 5, 13, 13, 13, 13 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum1 );                    \
  b3 = _mm512_fmadd_pd ( u1, t1, b3 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 6, 6, 6, 6, 10, 10, 10, 10 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b3 = _mm512_fmadd_pd ( u1, c2, b3 );                                  \
  indexes = _mm512_setr_epi64( 7, 7, 7, 7, 11, 11, 11, 11 );            \
  u1 = _mm512_permutex2var_pd( tu2, indexes, tum2 );                    \
  b3 = _mm512_fmadd_pd ( u1, t2, b3 );                                  \
                                                                        \
  indexes = _mm512_setr_epi64( 0, 0, 0, 0, 8, 8, 8, 8 );                \
  u1 = _mm512_permutex2var_pd( tu3, indexes, tum3 );                    \
  b3 = _mm512_fmadd_pd ( u1, c3, b3 );                                  \
  indexes = _mm512_setr_epi64( 1, 1, 1, 1, 9, 9, 9, 9 );                \
  u1 = _mm512_permutex2var_pd( tu3, indexes, tum3 );                    \
  b3 = _mm512_fmadd_pd ( u1, t3, b3 );                                  \
}






/* Combine spinor entries into 2 weyl vectors
   stored in high and low entries of a spinor
 */
#define   _avx512_to_weyl_1( w, b ){                                 \
  __m512i indexes;                                                   \
  __m512d _t;                                                        \
  indexes = _mm512_setr_epi64( 4, 5, 6, 7, 0, 1, 2, 3 );             \
  _t = _mm512_permutexvar_pd( indexes, (b) );                        \
  w = _mm512_maskz_add_pd( 0b00001111, _t, (b) );                    \
  w = _mm512_mask_sub_pd( w, 0b11110000, _t, (b) );                  \
}

#define   _avx512_to_weyl_2( w, b ){                                 \
  __m512i indexes;                                                   \
  __m512d _t;                                                        \
  indexes = _mm512_setr_epi64( 4, 5, 6, 7, 0, 1, 2, 3 );             \
  _t = _mm512_permutexvar_pd( indexes, (b) );                        \
  _t = _mm512_mask_add_pd( _t, 0b00001111, _t, (b) );                \
  _t = _mm512_mask_sub_pd( _t, 0b11110000, (b), _t );                \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 7, 6, 5, 4 );             \
  _t = _mm512_permutexvar_pd( indexes, _t );                         \
  w = _mm512_mask_add_pd( w, 0b10101111, w, _t );                    \
  w = _mm512_mask_sub_pd( w, 0b01010000, w, _t );                    \
}

#define   _avx512_to_weyl_3( w, b ){                                \
  __m512i indexes;                                                  \
  __m512d _t;                                                       \
  indexes = _mm512_setr_epi64( 4, 5, 6, 7, 0, 1, 2, 3 );            \
  _t = _mm512_permutexvar_pd( indexes, (b) );                       \
  _t = _mm512_mask_add_pd( _t, 0b00001111, _t, (b) );               \
  _t = _mm512_mask_sub_pd( _t, 0b11110000, (b), _t );               \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 6, 7, 4, 5 );            \
  _t = _mm512_permutexvar_pd( indexes, _t );                        \
  w = _mm512_mask_add_pd( w, 0b00111111, w, _t );                   \
  w = _mm512_mask_sub_pd( w, 0b11000000, w, _t );                   \
}

#define   _avx512_to_weyl_4( w, b ){                                \
  __m512i indexes;                                                  \
  __m512d _t;                                                       \
  indexes = _mm512_setr_epi64( 4, 5, 6, 7, 0, 1, 2, 3 );            \
  _t = _mm512_permutexvar_pd( indexes, (b) );                       \
  _t = _mm512_mask_add_pd( _t, 0b00001111, _t, (b) );               \
  _t = _mm512_mask_sub_pd( _t, 0b11110000, (b), _t );               \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 5, 4, 7, 6 );            \
  _t = _mm512_permutexvar_pd( indexes, _t );                        \
  w = _mm512_mask_add_pd( w, 0b01101111, w, _t );                   \
  w = _mm512_mask_sub_pd( w, 0b10010000, w, _t );                   \
}




/* Create a full Dirac vector by adding and subtracting the indeces of
 * a weyl vector */
#define _avx512_expand_weyl( a, w ){                                \
  __m512i indexes;                                                  \
  __m512d _t;                                                       \
  indexes = _mm512_setr_epi64( 4, 5, 6, 7, 0, 1, 2, 3 );            \
  _t = _mm512_permutexvar_pd( indexes, (w) );                       \
  a = _mm512_maskz_add_pd(   0b00001111, _t, w );                   \
  a = _mm512_mask_sub_pd( a, 0b11110000, _t, w );                   \
}

#define _avx512_expand_weyl_2( a, w ){                              \
  __m512i indexes;                                                  \
  __m512d _t1, _t2;                                                 \
  indexes = _mm512_setr_epi64( 7, 6, 5, 4, 7, 6, 5, 4 );            \
  _t1 = _mm512_permutexvar_pd( indexes, (w) );                      \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 0, 1, 2, 3 );            \
  _t2 = _mm512_permutexvar_pd( indexes, (w) );                      \
  a = _mm512_maskz_add_pd(   0b01011010, _t2, _t1 );                \
  a = _mm512_mask_sub_pd( a, 0b10100101, _t2, _t1 );                \
}

#define _avx512_expand_weyl_3( a, w ){                              \
  __m512i indexes;                                                  \
  __m512d _t1, _t2;                                                 \
  indexes = _mm512_setr_epi64( 6, 7, 4, 5, 6, 7, 4, 5 );            \
  _t1 = _mm512_permutexvar_pd( indexes, (w) );                      \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 0, 1, 2, 3 );            \
  _t2 = _mm512_permutexvar_pd( indexes, (w) );                      \
  a = _mm512_maskz_add_pd(   0b11000011, _t2, _t1 );                \
  a = _mm512_mask_sub_pd( a, 0b00111100, _t2, _t1 );                \
}

#define _avx512_expand_weyl_4( a, w ){                              \
  __m512i indexes;                                                  \
  __m512d _t1, _t2;                                                 \
  indexes = _mm512_setr_epi64( 5, 4, 7, 6, 5, 4, 7, 6 );            \
  _t1 = _mm512_permutexvar_pd( indexes, (w) );                      \
  indexes = _mm512_setr_epi64( 0, 1, 2, 3, 0, 1, 2, 3 );            \
  _t2 = _mm512_permutexvar_pd( indexes, (w) );                      \
  a = _mm512_maskz_add_pd(   0b10010110, _t2, _t1 );                \
  a = _mm512_mask_sub_pd( a, 0b01101001, _t2, _t1 );                \
}





/* Load four complex numbers. */
#define _avx512_load_4_d( v,  c1, c2, c3, c4 )                      \
{                                                                   \
  __m128d t128l, t128u;                                             \
  __m256d t256l ,t256u;                                             \
  t128l = _mm_loadu_pd( &(c1).re );                                 \
  t128u = _mm_loadu_pd( &(c2).re );                                 \
  t256l = _mm256_castpd128_pd256( t128l );                          \
  t256l = _mm256_insertf128_pd( t256l, t128u, 1 );                  \
  t128l = _mm_loadu_pd( &(c3).re );                                 \
  t128u = _mm_loadu_pd( &(c4).re );                                 \
  t256u = _mm256_castpd128_pd256( t128l );                          \
  t256u = _mm256_insertf128_pd( t256u, t128u, 1 );                  \
  v = _mm512_castpd256_pd512( t256l );                              \
  v = _mm512_insertf64x4( v, t256u, 1 );                            \
}

/* Store four complex numbers */
#define _avx512_store_4_d( r,  v1, v2, v3, v4 )                     \
{                                                                   \
  __m256d t256;                                                     \
  __m128d t128;                                                     \
  t256 = _mm512_extractf64x4_pd( r, 1 );                            \
  t128 = _mm256_extractf128_pd( t256, 1 );                          \
  _mm_storeu_pd( &(v4).re, t128  );                                 \
  t128 = _mm256_castpd256_pd128( t256 );                            \
  _mm_storeu_pd( &(v3).re, t128  );                                 \
  t256 = _mm512_castpd512_pd256( r );                               \
  t128 = _mm256_extractf128_pd( t256, 1 );                          \
  _mm_storeu_pd( &(v2).re, t128  );                                 \
  t128 = _mm256_castpd256_pd128( t256 );                            \
  _mm_storeu_pd( &(v1).re, t128  );                                 \
}




/* Adding a vectors to a spinors */

/* Load half spinors and organize for adding */
#define _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm )   \
  s1 = _mm512_loadu_pd( (sm) );                                \
  s2 = _mm512_loadu_pd( (sm)+8 );                              \
  s3 = _mm512_loadu_pd( (sp) );                                \
  s4 = _mm512_loadu_pd( (sp)+8 );                              \
                                                               \
  idx = _mm512_setr_epi64( 0,1,2,3, 8,9,10,11 );               \
  t1 = _mm512_permutex2var_pd( s1, idx, s3 );                  \
  idx = _mm512_setr_epi64( 4,5,6,7, 12,13,14,15 );             \
  t2 = _mm512_permutex2var_pd( s1, idx, s3 );                  \
  idx = _mm512_setr_epi64( 0,1,2,3, 8,9,10,11 );               \
  t3 = _mm512_permutex2var_pd( s2, idx, s4 );                  \

/* reorganize weyl vectors for adding */
#define _avx512_reorganize_a_ud_d( b1,b2,b3, a1,a2,a3 )        \
  idx = _mm512_setr_epi64( 0,1,8,9, 4,5,12,13 );               \
  a1 = _mm512_permutex2var_pd( b1, idx, b2 );                  \
  idx = _mm512_setr_epi64( 0,1,10,11, 4,5,14,15 );             \
  a2 = _mm512_permutex2var_pd( b3, idx, b1 );                  \
  idx = _mm512_setr_epi64( 2,3,10,11, 6,7,14,15 );             \
  a3 = _mm512_permutex2var_pd( b2, idx, b3 );                  \

/* store after adding */
#define _avx512_write_a_ud_d( a1, a2, a3, sp, sm ){            \
  __m256d l;                                                   \
  l = _mm512_castpd512_pd256( a1 );                            \
  _mm256_storeu_pd( (sm), l );                                  \
  l = _mm512_castpd512_pd256( a2 );                            \
  _mm256_storeu_pd( (sm)+4, l );                                \
  l = _mm512_castpd512_pd256( a3 );                            \
  _mm256_storeu_pd( (sm)+8, l );                                \
                                                               \
  l = _mm512_extractf64x4_pd( a1, 1 );                         \
  _mm256_storeu_pd( (sp), l );                                  \
  l = _mm512_extractf64x4_pd( a2, 1 );                         \
  _mm256_storeu_pd( (sp)+4, l );                                \
  l = _mm512_extractf64x4_pd( a3, 1 );                         \
  _mm256_storeu_pd( (sp)+8, l );                                \
}


#define _avx512_add_to_spinors( b1, b2, b3, sp, sm )           \
{                                                              \
  __m512d s1,s2,s3,s4, a1,a2,a3, t1,t2,t3;                     \
  __m512i idx;                                                 \
                                                               \
  _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm );        \
  _avx512_reorganize_a_ud_d( b1,b2,b3, a1,a2,a3 );             \
                                                               \
  t1 = _mm512_add_pd( a1, t1 );                                \
  t2 = _mm512_add_pd( a2, t2 );                                \
  t3 = _mm512_add_pd( a3, t3 );                                \
                                                               \
  _avx512_write_a_ud_d( t1, t2, t3, sp, sm );                  \
}

#define _avx512_add_to_spinors_2( b1, b2, b3, sp, sm )         \
{                                                              \
  __m512d s1,s2,s3,s4, a1,a2,a3, t1,t2,t3;                     \
  __m512i idx;                                                 \
                                                               \
  _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm );        \
  _avx512_reorganize_a_ud_d( b1,b2,b3, a1,a2,a3 );             \
                                                               \
  t1 = _mm512_mask_add_pd( t1, 0b00001111, t1, a1 );           \
  t1 = _mm512_mask_sub_pd( t1, 0b11110000, t1, a1 );           \
  t2 = _mm512_mask_add_pd( t2, 0b00001111, t2, a2 );           \
  t2 = _mm512_mask_sub_pd( t2, 0b11110000, t2, a2 );           \
  t3 = _mm512_mask_add_pd( t3, 0b00001111, t3, a3 );           \
  t3 = _mm512_mask_sub_pd( t3, 0b11110000, t3, a3 );           \
                                                               \
  _avx512_write_a_ud_d( t1, t2, t3, sp, sm );                  \
}

#define _avx512_add_to_spinors_3( b1, b2, b3, sp, sm )         \
{                                                              \
  __m512d s1,s2,s3,s4, a1,a2,a3, t1,t2,t3;                     \
  __m512i idx;                                                 \
                                                               \
  _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm );        \
  idx = _mm512_setr_epi64( 3,2,11,10, 7,6,15,14 );             \
  a1 = _mm512_permutex2var_pd( b1, idx, b2 );                  \
  idx = _mm512_setr_epi64( 3,2,9,8, 7,6,13,12 );               \
  a2 = _mm512_permutex2var_pd( b3, idx, b1 );                  \
  idx = _mm512_setr_epi64( 1,0,9,8, 5,4,13,12 );               \
  a3 = _mm512_permutex2var_pd( b2, idx, b3 );                  \
                                                               \
  t1 = _mm512_mask_add_pd( t1, 0b10100101, t1, a1 );           \
  t1 = _mm512_mask_sub_pd( t1, 0b01011010, t1, a1 );           \
  t2 = _mm512_mask_add_pd( t2, 0b10100101, t2, a2 );           \
  t2 = _mm512_mask_sub_pd( t2, 0b01011010, t2, a2 );           \
  t3 = _mm512_mask_add_pd( t3, 0b10100101, t3, a3 );           \
  t3 = _mm512_mask_sub_pd( t3, 0b01011010, t3, a3 );           \
                                                               \
  _avx512_write_a_ud_d( t1, t2, t3, sp, sm );                  \
}

#define _avx512_add_to_spinors_4( b1, b2, b3, sp, sm )         \
{                                                              \
  __m512d s1,s2,s3,s4, a1,a2,a3, t1,t2,t3;                     \
  __m512i idx;                                                 \
  idx = _mm512_setr_epi64( 2,3,0,1, 6,7,4,5 );                 \
  b1 = _mm512_permutexvar_pd( idx, b1 );                       \
  b2 = _mm512_permutexvar_pd( idx, b2 );                       \
  b3 = _mm512_permutexvar_pd( idx, b3 );                       \
                                                               \
  _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm );        \
  idx = _mm512_setr_epi64( 0,1,8,9, 4,5,12,13 );               \
  a1 = _mm512_permutex2var_pd( b1, idx, b2 );                  \
  idx = _mm512_setr_epi64( 0,1,10,11, 4,5,14,15 );             \
  a2 = _mm512_permutex2var_pd( b3, idx, b1 );                  \
  idx = _mm512_setr_epi64( 2,3,10,11, 6,7,14,15 );             \
  a3 = _mm512_permutex2var_pd( b2, idx, b3 );                  \
                                                               \
  t1 = _mm512_mask_add_pd( t1, 0b11110000, t1, a1 );           \
  t1 = _mm512_mask_sub_pd( t1, 0b00001111, t1, a1 );           \
  t2 = _mm512_mask_add_pd( t2, 0b00111100, t2, a2 );           \
  t2 = _mm512_mask_sub_pd( t2, 0b11000011, t2, a2 );           \
  t3 = _mm512_mask_add_pd( t3, 0b00001111, t3, a3 );           \
  t3 = _mm512_mask_sub_pd( t3, 0b11110000, t3, a3 );           \
                                                               \
  _avx512_write_a_ud_d( t1, t2, t3, sp, sm );                  \
}


#define _avx512_add_to_spinors_5( b1, b2, b3, sp, sm )         \
{                                                              \
  __m512d s1,s2,s3,s4, a1,a2,a3, t1,t2,t3;                     \
  __m512i idx;                                                 \
                                                               \
  _avx512_load_s_ud_d( s1,s2,s3,s4, t1,t2,t3, sp, sm );        \
  idx = _mm512_setr_epi64( 1,0,9,8, 5,4,13,12 );               \
  a1 = _mm512_permutex2var_pd( b1, idx, b2 );                  \
  idx = _mm512_setr_epi64( 1,0,11,10, 5,4,15,14 );             \
  a2 = _mm512_permutex2var_pd( b3, idx, b1 );                  \
  idx = _mm512_setr_epi64( 3,2,11,10, 7,6,15,14 );             \
  a3 = _mm512_permutex2var_pd( b2, idx, b3 );                  \
                                                               \
  t1 = _mm512_mask_add_pd( t1, 0b10100101, t1, a1 );           \
  t1 = _mm512_mask_sub_pd( t1, 0b01011010, t1, a1 );           \
  t2 = _mm512_mask_add_pd( t2, 0b01101001, t2, a2 );           \
  t2 = _mm512_mask_sub_pd( t2, 0b10010110, t2, a2 );           \
  t3 = _mm512_mask_add_pd( t3, 0b01011010, t3, a3 );           \
  t3 = _mm512_mask_sub_pd( t3, 0b10100101, t3, a3 );           \
                                                               \
  _avx512_write_a_ud_d( t1, t2, t3, sp, sm );                  \
}





#endif
