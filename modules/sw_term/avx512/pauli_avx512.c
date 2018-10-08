/*******************************************************************************
*
* File pauli_avx512.c
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementations of the clover term multiplication in
* single precision.
*
* See ../pauli_avx512.c for more information and alternative
* implementations.
*
*******************************************************************************/

#ifdef AVX512

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "sw_term.h"
typedef union
{
  spinor s;
  weyl w[2];
} spin_t;

#include "avx512.h"


void mul_pauli2(float mu, pauli *m, spinor *source, spinor *res )
{
  spin_t *ps, *pr;
  float const *u, *u2;
  __m512i idx;
  __m512 tr1,tr2,tr3;
  __m512 ts1, ts2, ts3, tsi1, tsi2, tsi3, u512;
  __m512 tu11, tu12, tu13, tu21, tu22, tu23;
  __m512 tu1, tu2, tu3, tu4;
  __m512 umu;
  __m256 t256;
  __m128 t128a, t128b, tmu;

  ps = (spin_t *)(source);
  pr = (spin_t *)(res);

  u = (*m).u;
  u2 = (m+1)->u;

  weyl * s = (*ps).w;
  weyl * r = (*pr).w;
  weyl * s2 = (*ps).w+1;
  weyl * r2 = (*pr).w+1;

  s += 4;
  _prefetch_spinor(s);
  s -= 4;

  tr1 = _mm512_loadu_ps( &(*s).c1.c1.re );
  tr2 = _mm512_castps256_ps512( _mm256_loadu_ps( &(*s).c1.c1.re+16 ) );
  idx = _mm512_setr_epi32( 0,1,2,3,6,7,8,9, 12,13,14,15,18,19,20,21 );
  ts1 = _mm512_permutex2var_ps( tr1, idx, tr2 );
  idx = _mm512_setr_epi32( 2,3,4,5,8,9,10,11, 14,15,16,17,20,21,22,23 );
  ts2 = _mm512_permutex2var_ps( tr1, idx, tr2 );
  idx = _mm512_setr_epi32( 4,5,0,1,10,11,6,7, 16,17,12,13,22,23,18,19 );
  ts3 = _mm512_permutex2var_ps( tr1, idx, tr2 );

  tu11 = _mm512_loadu_ps( u );
  tu12 = _mm512_loadu_ps( u+16 );
  tu13 = _mm512_loadu_ps( u+32 );
  tu21 = _mm512_loadu_ps( u2 );
  tu22 = _mm512_loadu_ps( u2+16 );
  tu23 = _mm512_loadu_ps( u2+32 );


  tsi1 = _mm512_permute_ps ( ts1, 0b10110001 );
  tsi2 = _mm512_permute_ps ( ts2, 0b10110001 );
  tsi3 = _mm512_permute_ps ( ts3, 0b10110001 );


  tmu = _mm_load_ps1( &mu );
  umu = _mm512_broadcastss_ps( tmu );


  idx = _mm512_setr_epi32( 0,1,10,11,4,5,2,3, 16,17,26,27,20,21,18,19 );
  tu1 = _mm512_permutex2var_ps( tu11, idx, tu21 );
  idx = _mm512_setr_epi32( 4,5,0,1,12,13,6,7, 20,21,16,17,28,29,22,23);
  tu2 = _mm512_permutex2var_ps( tu12, idx, tu22 );

  idx = _mm512_setr_epi32( 0,0,1,1,2,3,16+0,16+1, 8,8,9,9,10,11,16+8,16+9 );
  tu4 = _mm512_permutex2var_ps( tu1, idx, tu2 );
  u512 = _mm512_permute_ps( tu4, 0b10100000 );
  tr1 = _mm512_mul_ps( ts1, u512 );

  idx = _mm512_setr_epi32( 0,1,2,3,16+5,16+5,16+7,16+7,
                           8,9,10,11,16+13,16+13,16+15,16+15 );
  u512 = _mm512_permutex2var_ps( umu, idx, tu4 );
  u512 = _mm512_mul_ps( u512, tsi1 );
  tr1 = _mm512_mask_add_ps( tr1, 0b1010010110101010, tr1, u512 );
  tr1 = _mm512_mask_sub_ps( tr1, 0b0101101001010101, tr1, u512 );


  idx = _mm512_setr_epi32( 0,0,4,4,16+4,16+4,16+5,16+5,
                           8,8,12,12,16+12,16+12,16+13,16+13 );
  u512 = _mm512_permutex2var_ps( tu2, idx, tu1 );
  tr2 = _mm512_mul_ps( ts2, u512 );

  idx = _mm512_setr_epi32( 1,1,5,5,16+4,16+5,16+6,16+7,
                           9,9,13,13,16+12,16+13,16+14,16+15 );
  u512 = _mm512_permutex2var_ps( tu2, idx, umu );
  u512 = _mm512_mul_ps( u512, tsi2 );
  tr2 = _mm512_mask_add_ps( tr2, 0b0101010110100101, tr2, u512 );
  tr2 = _mm512_mask_sub_ps( tr2, 0b1010101001011010, tr2, u512 );


  idx = _mm512_setr_epi32( 6,6,2,3,16+4,16+5,7,7,
                           14,14,10,11,16+12,16+13,15,15 );
  tu4 = _mm512_permutex2var_ps( tu1, idx, tu2 );
  u512 = _mm512_permute_ps( tu4, 0b10100000 );
  tr3 = _mm512_mul_ps( ts3, u512 );

  idx = _mm512_setr_epi32( 0,1,16+3,16+3,16+5,16+5,6,7,
                           8,9,16+11,16+11,16+13,16+13,14,15 );
  u512 = _mm512_permutex2var_ps( umu, idx, tu4 );
  u512 = _mm512_mul_ps( u512, tsi3 );
  tr3 = _mm512_mask_add_ps( tr3, 0b0110010110100110, tr3, u512 );
  tr3 = _mm512_mask_sub_ps( tr3, 0b1001101001011001, tr3, u512 );



  idx = _mm512_setr_epi32( 8,9,6,7,14,15,12,13,
                           24,25,22,23,30,31,28,29 );
  tu1 = _mm512_permutex2var_ps( tu11, idx, tu21 );

  u512 = _mm512_shuffle_ps( tu1, tu2, 0b10101010 );
  tr1 = _mm512_fmadd_ps( ts2, u512, tr1 );

  u512 = _mm512_shuffle_ps( tu1, tu2, 0b11111111 );
  u512 = _mm512_mul_ps( u512, tsi2 );
  tr1 = _mm512_mask_add_ps( tr1, 0b1010101010101010, tr1, u512 );
  tr1 = _mm512_mask_sub_ps( tr1, 0b0101010101010101, tr1, u512 );

  m += 4;
  _prefetch_pauli_dble(m);
  m -= 4;

  idx = _mm512_setr_epi32( 0,1,2,3,0,1,2,3, 16,17,18,19,16,17,18,19 );
  tu13 = _mm512_permutex2var_ps( tu13, idx, tu23 );
  idx = _mm512_setr_epi32( 6,7,2,3,8,9,14,15, 22,23,18,19,24,25,30,31 );
  tu2 = _mm512_permutex2var_ps( tu12, idx, tu22 );

  idx = _mm512_setr_epi32( 6,7,16+0,16+1,16+6,16+7,0,0,
                           14,15,16+8,16+9,16+14,16+15,0,0 );
  tu3 = _mm512_permutex2var_ps( tu1, idx, tu2 );
  tu4 = _mm512_permute_ps( tu3, 0b10100000 );
  u512 = _mm512_mask_shuffle_ps( tu4, 0b1111000011110000, tu4, tu13, 0b10100100 );
  tr2 = _mm512_fmadd_ps( ts1, u512, tr2 );

  tu4 = _mm512_permute_ps( tu3, 0b11110101 );
  u512 = _mm512_mask_shuffle_ps( tu4, 0b1111000011110000, tu4, tu13, 0b11110100 );
  u512 = _mm512_mul_ps( u512, tsi1 );
  tr2 = _mm512_mask_add_ps( tr2, 0b0101010101010101, tr2, u512 );
  tr2 = _mm512_mask_sub_ps( tr2, 0b1010101010101010, tr2, u512 );

  tu3 = _mm512_mask_shuffle_ps( tu2, 0b0000111100001111, tu1, tu2, 0b11100100 );
  u512 = _mm512_permute_ps( tu3, 0b10100000 );
  tr3 = _mm512_fmadd_ps( ts1, u512, tr3 );

  u512 = _mm512_permute_ps( tu3, 0b11110101 );
  u512 = _mm512_mul_ps( u512, tsi1 );
  tr3 = _mm512_mask_add_ps( tr3, 0b1010010110100101, tr3, u512 );
  tr3 = _mm512_mask_sub_ps( tr3, 0b0101101001011010, tr3, u512 );

  idx = _mm512_setr_epi32( 0,1,2,3, 4,5,16+2,16+3, 8,9,10,11,12,13,16+10,16+11 );
  tu3 = _mm512_permutex2var_ps( tu1, idx, tu2 );
  u512 = _mm512_permute_ps( tu3, 0b10100000 );
  tr1 = _mm512_fmadd_ps( ts3, u512, tr1 );

  u512 = _mm512_permute_ps( tu3, 0b11110101 );
  u512 = _mm512_mul_ps( u512, tsi3 );
  tr1 = _mm512_mask_add_ps( tr1, 0b1010011010100110, tr1, u512 );
  tr1 = _mm512_mask_sub_ps( tr1, 0b0101100101011001, tr1, u512 );


  idx = _mm512_setr_epi32( 0,1,8,9,10,11,-1,-1, 16,17,24,25,26,27,-1,-1 );
  tu2 = _mm512_permutex2var_ps( tu12, idx, tu22 );

  idx = _mm512_setr_epi32( 4,5,16+4,16+5,0,0,0,0, 12,13,16+12,16+13,0,0,0,0 );
  tu3 = _mm512_permutex2var_ps( tu2, idx, tu1 );
  u512 = _mm512_permute_ps( tu3, 0b10100000 );
  u512 = _mm512_mask_permute_ps( u512, 0b1111000011110000, tu13, 0b00001010 );
  tr2 = _mm512_fmadd_ps( ts3, u512, tr2 );

  u512 = _mm512_permute_ps( tu3, 0b11110101 );
  u512 = _mm512_mask_permute_ps( u512, 0b1111000011110000, tu13, 0b01011111 );
  u512 = _mm512_mul_ps( u512, tsi3 );
  tr2 = _mm512_mask_add_ps( tr2, 0b0110010101100101, tr2, u512 );
  tr2 = _mm512_mask_sub_ps( tr2, 0b1001101010011010, tr2, u512 );

  tu3 = _mm512_mask_shuffle_ps( tu2, 0b1111000011110000, tu2, tu13, 0b01000100 );
  u512 = _mm512_permute_ps( tu3, 0b10100000 );
  tr3 = _mm512_fmadd_ps( ts2, u512, tr3 );

  u512 = _mm512_permute_ps( tu3, 0b11110101 );
  u512 = _mm512_mul_ps( u512, tsi2 );
  tr3 = _mm512_mask_add_ps( tr3, 0b1010010110100101, tr3, u512 );
  tr3 = _mm512_mask_sub_ps( tr3, 0b0101101001011010, tr3, u512 );



  idx = _mm512_setr_epi32( 0,1,2,3, 16,17,18,19, 8,9,10,11, 24,25,26,27 );
  ts1 = _mm512_permutex2var_ps( tr1, idx, tr3 );
  idx = _mm512_setr_epi32( 4,5,6,7, 20,21,22,23, 12,13,14,15, 28,29,30,31 );
  ts2 = _mm512_permutex2var_ps( tr1, idx, tr3 );
  ts3 = _mm512_add_ps( ts1, ts2 );

  t256 = _mm512_castps512_ps256( ts3 );
  _mm256_storeu_ps( &(*r).c1.c1.re, t256 );

  t128a = _mm512_castps512_ps128( tr2 );
  t128b = _mm512_extractf32x4_ps( tr2, 1 );
  t128b = _mm_add_ps( t128a, t128b );
  _mm_storeu_ps( &(*r).c2.c2.re, t128b );

  t256 = _mm256_castpd_ps( _mm512_extractf64x4_pd( _mm512_castps_pd(ts3), 1 ) );
  _mm256_storeu_ps( &(*r2).c1.c1.re, t256 );

  t128a = _mm512_extractf32x4_ps( tr2, 2 );
  t128b = _mm512_extractf32x4_ps( tr2, 3 );
  t128b = _mm_add_ps( t128a, t128b );
  _mm_storeu_ps( &(*r2).c2.c2.re, t128b );
}

#endif