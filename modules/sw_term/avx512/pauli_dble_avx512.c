/*******************************************************************************
*
* File pauli_avx512.c
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementations of the clover term multiplication in
* double precision.
*
* See ../pauli_dble_avx512.c for more information and alternative
* implementations.
*
*******************************************************************************/

#ifdef AVX512

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "su3.h"
#include "linalg.h"
#include "sw_term.h"
#define DELTA 1.0e-04

typedef union
{
  spinor_dble s;
  weyl_dble w[2];
  complex_dble c[12];
} spin_t;

#include "avx512.h"


void mul_pauli2_dble(double mu, pauli_dble *m, weyl_dble *s, weyl_dble *r)
{
  double const  *u = m->u, *u2 = (m+1)->u;

  __m512d r1,r2,r3;
  __m512d s1,s2,s3,s4,s5,s6, si1,si2,si3,si4,si5,si6, s512, u512;
  __m512d t1,t2,t3,t4,t5,t6;
  __m512d tu11, tu12, tu13, tu14, tu15, tu21, tu22, tu23, tu24, tu25;
  __m512d tu1, tu2, tu3, tu4, tu5, tu6;
  __m512d umu;
  __m512i idx;

  t1 = _mm512_loadu_pd( &(*s).c1.c1.re );
  t2 = _mm512_loadu_pd( &(*s).c1.c1.re + 8 );
  t3 = _mm512_loadu_pd( &(*s).c1.c1.re + 16 );

  tu11 = _mm512_loadu_pd( u );
  tu12 = _mm512_loadu_pd( u+8 );
  tu13 = _mm512_loadu_pd( u+16 );
  tu14 = _mm512_loadu_pd( u+24 );
  tu15 = _mm512_loadu_pd( u+32 );
  tu21 = _mm512_loadu_pd( u2 );
  tu22 = _mm512_loadu_pd( u2+8 );
  tu23 = _mm512_loadu_pd( u2+16 );
  tu24 = _mm512_loadu_pd( u2+24 );
  tu25 = _mm512_loadu_pd( u2+32 );

  umu = _mm512_loadu_pd( &mu );


  idx = _mm512_setr_epi64( 0,1,12,13,0,1,12,13 );
  s1 = _mm512_permutex2var_pd( t1, idx, t2 );
  idx = _mm512_setr_epi64( 2,3,14,15,2,3,14,15 );
  s2 = _mm512_permutex2var_pd( t1, idx, t2 );
  idx = _mm512_setr_epi64( 4,5,8,9,4,5,8,9 );
  s3 = _mm512_permutex2var_pd( t1, idx, t3 );
  idx = _mm512_setr_epi64( 0,1,12,13,0,1,12,13 );
  s4 = _mm512_permutex2var_pd( t2, idx, t3 );
  idx = _mm512_setr_epi64( 6,7,10,11,6,7,10,11 );
  s5 = _mm512_permutex2var_pd( t1, idx, t3 );
  idx = _mm512_setr_epi64( 2,3,14,15,2,3,14,15 );
  s6 = _mm512_permutex2var_pd( t2, idx, t3 );

  si1 = _mm512_permute_pd ( s1, 0b01010101 );
  si2 = _mm512_permute_pd ( s2, 0b01010101 );
  si3 = _mm512_permute_pd ( s3, 0b01010101 );
  si4 = _mm512_permute_pd ( s4, 0b01010101 );
  si5 = _mm512_permute_pd ( s5, 0b01010101 );
  si6 = _mm512_permute_pd ( s6, 0b01010101 );


  idx = _mm512_setr_epi64( 0,1,8,9,6,7,14,15 );
  tu1 = _mm512_permutex2var_pd( tu11, idx, tu21 );
  idx = _mm512_setr_epi64( 4,5,12,13,2,3,10,11 );
  tu2 = _mm512_permutex2var_pd( tu12, idx, tu22 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu2 );
  r1 = _mm512_mul_pd( u512, s1 );
  idx = _mm512_setr_epi64( 0,0,0,0,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( umu, idx, tu2 );
  u512 = _mm512_mul_pd( u512, si1 );
  r1 = _mm512_mask_add_pd( r1, 0b01010110, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b10101001, r1, u512 );

  idx = _mm512_setr_epi64( 2,3,10,11,4,5,12,13 );
  tu6 = _mm512_permutex2var_pd( tu11, idx, tu21 );

  idx = _mm512_setr_epi64( 4,4,6,6,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu6 );
  r1 = _mm512_fmadd_pd( u512, s5, r1 );
  idx = _mm512_setr_epi64( 5,5,7,7,8,8,8,8 );
  u512 = _mm512_permutex2var_pd( tu2, idx, umu );
  u512 = _mm512_mul_pd( u512, si5 );
  r1 = _mm512_mask_add_pd( r1, 0b01101010, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b10010101, r1, u512 );


  idx = _mm512_setr_epi64( 2,3,8+2,8+3,4,5,8+4,8+5 );
  tu3 = _mm512_permutex2var_pd( tu13, idx, tu23 );

  idx = _mm512_setr_epi64( 1,1,3,3,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu3 );
  r2 = _mm512_mul_pd( u512, s2 );
  idx = _mm512_setr_epi64( 0,0,0,0,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( umu, idx, tu3 );
  u512 = _mm512_mul_pd( u512, si2 );
  r2 = _mm512_mask_add_pd( r2, 0b01010110, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b10101001, r2, u512 );

  idx = _mm512_setr_epi64( 4,4,6,6,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu6 );
  r2 = _mm512_fmadd_pd( u512, s4, r2 );
  idx = _mm512_setr_epi64( 5,5,7,7,8,8,8,8 );
  u512 = _mm512_permutex2var_pd( tu3, idx, umu );
  u512 = _mm512_mul_pd( u512, si4 );
  r2 = _mm512_mask_add_pd( r2, 0b01101010, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b10010101, r2, u512 );


  idx = _mm512_setr_epi64( 4,5,8+4,8+5,6,7,8+6,8+7 );
  tu4 = _mm512_permutex2var_pd( tu14, idx, tu24 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+0,8+0,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu6, idx, tu4 );
  r3 = _mm512_mul_pd( u512, s3 );
  idx = _mm512_setr_epi64( 0,0,0,0,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( umu, idx, tu4 );
  u512 = _mm512_mul_pd( u512, si3 );
  r3 = _mm512_mask_add_pd( r3, 0b01010110, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10101001, r3, u512 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu4, idx, tu6 );
  r3 = _mm512_fmadd_pd( u512, s6, r3 );
  idx = _mm512_setr_epi64( 1,1,3,3,8,8,8,8 );
  u512 = _mm512_permutex2var_pd( tu4, idx, umu );
  u512 = _mm512_mul_pd( u512, si6 );
  r3 = _mm512_mask_add_pd( r3, 0b01101010, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10010101, r3, u512 );


  idx = _mm512_setr_epi64( 4,4,6,6,8,8,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu2 );
  r2 = _mm512_fmadd_pd( u512, s1, r2 );
  idx = _mm512_setr_epi64( 5,5,7,7,9,9,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu2 );
  u512 = _mm512_mul_pd( u512, si1 );
  r2 = _mm512_mask_add_pd( r2, 0b01010101, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b10101010, r2, u512 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu4 );
  r1 = _mm512_fmadd_pd( u512, s4, r1 );
  idx = _mm512_setr_epi64( 1,1,3,3,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu4 );
  u512 = _mm512_mul_pd( u512, si4 );
  r1 = _mm512_mask_add_pd( r1, 0b10101010, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b01010101, r1, u512 );


  idx = _mm512_setr_epi64( 4,4,6,6,8+0,8+0,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu3 );
  r1 = _mm512_fmadd_pd( u512, s2, r1 );
  idx = _mm512_setr_epi64( 5,5,7,7,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu1, idx, tu3 );
  u512 = _mm512_mul_pd( u512, si2 );
  r1 = _mm512_mask_add_pd( r1, 0b01011010, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b10100101, r1, u512 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu4 );
  r2 = _mm512_fmadd_pd( u512, s5, r2 );
  idx = _mm512_setr_epi64( 1,1,3,3,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu4 );
  u512 = _mm512_mul_pd( u512, si5 );
  r2 = _mm512_mask_add_pd( r2, 0b01011010, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b10100101, r2, u512 );

  idx = _mm512_setr_epi64( 0,0,8,8,6,6,14,14 );
  u512 = _mm512_permutex2var_pd( tu12, idx, tu22 );
  r3 = _mm512_fmadd_pd( u512, s1, r3 );
  idx = _mm512_setr_epi64( 1,1,9,9,7,7,15,15 );
  u512 = _mm512_permutex2var_pd( tu12, idx, tu22 );
  u512 = _mm512_mul_pd( u512, si1 );
  r3 = _mm512_mask_add_pd( r3, 0b01010101, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10101010, r3, u512 );


  idx = _mm512_setr_epi64( 0,0,8,8,6,6,14,14 );
  u512 = _mm512_permutex2var_pd( tu13, idx, tu23 );
  r3 = _mm512_fmadd_pd( u512, s2, r3 );
  idx = _mm512_setr_epi64( 1,1,9,9,7,7,15,15 );
  u512 = _mm512_permutex2var_pd( tu13, idx, tu23 );
  u512 = _mm512_mul_pd( u512, si2 );
  r3 = _mm512_mask_add_pd( r3, 0b01010101, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10101010, r3, u512 );



  idx = _mm512_setr_epi64( 0,1,8+0,8+1,6,7,8+6,8+7 );
  tu2 = _mm512_permutex2var_pd( tu12, idx, tu22 );
  idx = _mm512_setr_epi64( 2,3,8+2,8+3,0,1,8+0,8+1 );
  tu4 = _mm512_permutex2var_pd( tu14, idx, tu24 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu4 );
  r1 = _mm512_fmadd_pd( u512, s3, r1 );
  idx = _mm512_setr_epi64( 1,1,3,3,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu4 );
  u512 = _mm512_mul_pd( u512, si3 );
  r1 = _mm512_mask_add_pd( r1, 0b01011010, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b10100101, r1, u512 );


  idx = _mm512_setr_epi64( 0,1,8+0,8+1,2,3,8+2,8+3 );
  tu5 = _mm512_permutex2var_pd( tu15, idx, tu25 );

  idx = _mm512_setr_epi64( 4,4,6,6,8+0,8+0,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu4, idx, tu5 );
  r3 = _mm512_fmadd_pd( u512, s5, r3 );
  idx = _mm512_setr_epi64( 5,5,7,7,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu4, idx, tu5 );
  u512 = _mm512_mul_pd( u512, si5 );
  r3 = _mm512_mask_add_pd( r3, 0b01011010, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10100101, r3, u512 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu4, idx, tu5 );
  r3 = _mm512_fmadd_pd( u512, s4, r3 );
  idx = _mm512_setr_epi64( 1,1,3,3,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu4, idx, tu5 );
  u512 = _mm512_mul_pd( u512, si4 );
  r3 = _mm512_mask_add_pd( r3, 0b01011010, r3, u512 );
  r3 = _mm512_mask_sub_pd( r3, 0b10100101, r3, u512 );

  idx = _mm512_setr_epi64( 4,4,6,6,8+0,8+0,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu5 );
  r1 = _mm512_fmadd_pd( u512, s6, r1 );
  idx = _mm512_setr_epi64( 5,5,7,7,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu2, idx, tu5 );
  u512 = _mm512_mul_pd( u512, si6 );
  r1 = _mm512_mask_add_pd( r1, 0b10101010, r1, u512 );
  r1 = _mm512_mask_sub_pd( r1, 0b01010101, r1, u512 );


  idx = _mm512_setr_epi64( 6,7,8+6,8+7,0,1,8+0,8+1 );
  tu3 = _mm512_permutex2var_pd( tu13, idx, tu23 );

  idx = _mm512_setr_epi64( 4,4,6,6,8+0,8+0,8+2,8+2 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu4 );
  r2 = _mm512_fmadd_pd( u512, s3, r2 );
  idx = _mm512_setr_epi64( 5,5,7,7,8+1,8+1,8+3,8+3 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu4 );
  u512 = _mm512_mul_pd( u512, si3 );
  r2 = _mm512_mask_add_pd( r2, 0b01011010, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b10100101, r2, u512 );

  idx = _mm512_setr_epi64( 0,0,2,2,8+4,8+4,8+6,8+6 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu5 );
  r2 = _mm512_fmadd_pd( u512, s6, r2 );
  idx = _mm512_setr_epi64( 1,1,3,3,8+5,8+5,8+7,8+7 );
  u512 = _mm512_permutex2var_pd( tu3, idx, tu5 );
  u512 = _mm512_mul_pd( u512, si6 );
  r2 = _mm512_mask_add_pd( r2, 0b10101010, r2, u512 );
  r2 = _mm512_mask_sub_pd( r2, 0b01010101, r2, u512 );

  idx = _mm512_setr_epi64( 0,1,8,9,2,3,10,11 );
  t1 = _mm512_permutex2var_pd( r1, idx, r2 );
  idx = _mm512_setr_epi64( 2,3,14,15,0,1,12,13 );
  t2 = _mm512_permutex2var_pd( r3, idx, r1 );
  idx = _mm512_setr_epi64( 4,5,12,13,6,7,14,15 );
  t3 = _mm512_permutex2var_pd( r2, idx, r3 );
  r1 = _mm512_mask_blend_pd( 0b11110000, t1, t2 );
  r2 = _mm512_mask_blend_pd( 0b11110000, t3, t1 );
  r3 = _mm512_mask_blend_pd( 0b11110000, t2, t3 );

  _mm512_storeu_pd( &r[0].c1.c1.re, r1 );
  _mm512_storeu_pd( &r[0].c2.c2.re, r2 );
  _mm512_storeu_pd( &r[1].c1.c3.re, r3 );
}


int fwd_house_avx512(double eps, complex_dble *aa, complex_dble *dd, double * rr )
{
  int i, j, k, ifail;
  double r1, r2, r3;
  complex_dble z;

  ifail = 0;

  for (k = 0; k < 5; k++) {
    r1 = aa[6 * k + k].re * aa[6 * k + k].re +
         aa[6 * k + k].im * aa[6 * k + k].im;
    r2 = sqrt(r1);

    for (j = (k + 1); j < 6; j++)
      r1 += (aa[6 * j + k].re * aa[6 * j + k].re +
             aa[6 * j + k].im * aa[6 * j + k].im);

    if (r1 >= eps)
      r1 = sqrt(r1);
    else {
      ifail = 1;
      r1 = 1.0;
    }

    if (r2 >= (DBL_EPSILON * r1)) {
      r3 = 1.0 / r2;
      z.re = r3 * aa[6 * k + k].re;
      z.im = r3 * aa[6 * k + k].im;
    } else {
      z.re = 1.0;
      z.im = 0.0;
    }

    aa[6 * k + k].re += r1 * z.re;
    aa[6 * k + k].im += r1 * z.im;

    r3 = 1.0 / (r1 * (r1 + r2));
    rr[k] = r3;
    dd[k].re = -(r1 + r2) * r3 * z.re;
    dd[k].im = (r1 + r2) * r3 * z.im;

    for (j = (k + 1); j < 6; j++) {
      complex_dble z, *ak, *aj;
      __m128d mz, t1, t2, t3;
      mz = _mm_setzero_pd();

      ak = aa + 6 * k + k;
      aj = aa + 6 * k + j;

      for (i = k; i < 6; i++) {
        t1 = _mm_loaddup_pd(&ak->re);
        t2 = _mm_loaddup_pd(&ak->im);
        t3 = _mm_load_pd(&aj->re);
        t2 = _mm_mul_pd( t2, t3 );
        t2 = _mm_permute_pd( t2, 0b01 );
        t1 = _mm_fmsubadd_pd( t1, t3, t2 );
        mz = _mm_add_pd( mz, t1 );
        ak += 6;
        aj += 6;
      }

      t1 = _mm_loaddup_pd(&r3);
      mz = _mm_mul_pd( mz, t1 );
      _mm_storeu_pd( &z.re, mz );

      ak = aa + 6 * k + k;
      aj = aa + 6 * k + j;
      for (i = k; i < 6; i++) {
        t1 = _mm_loaddup_pd(&ak->re);
        t2 = _mm_loaddup_pd(&ak->im);
        t3 = _mm_load_pd(&aj->re);
        t2 = _mm_mul_pd( mz, t2 );
        t2 = _mm_permute_pd( t2, 0b01 );
        t1 = _mm_fmaddsub_pd( mz,t1, t2 );
        t3 = _mm_sub_pd( t3, t1 );
        _mm_storeu_pd( &aj->re, t3 );
        ak += 6;
        aj += 6;
      }
    }
  }

  r1 = aa[35].re * aa[35].re + aa[35].im * aa[35].im;

  if (r1 >= eps)
    r1 = 1.0 / r1;
  else {
    ifail = 1;
    r1 = 1.0;
  }

  dd[5].re = r1 * aa[35].re;
  dd[5].im = -r1 * aa[35].im;

  return ifail;
}


void solv_sys_avx512( complex_dble *aa, complex_dble *dd )
{
  int i, j, k;
  complex_dble z;
  __m128d mz, t1, t2, t3;

  for (k = 5; k > 0; k--) {
    for (i = (k - 1); i >= 0; i--) {
      t1 = _mm_loaddup_pd(&aa[6 * i + k].re);
      t2 = _mm_loaddup_pd(&aa[6 * i + k].im);
      t3 = _mm_load_pd(&dd[k].re);
      t2 = _mm_mul_pd( t2, t3 );
      t2 = _mm_permute_pd( t2, 0b01 );
      mz = _mm_fmaddsub_pd( t1, t3, t2 );

      for (j = (k - 1); j > i; j--) {
        t1 = _mm_loaddup_pd(&aa[6 * i + j].re);
        t2 = _mm_loaddup_pd(&aa[6 * i + j].im);
        t3 = _mm_load_pd(&aa[6 * j + k].re);
        t2 = _mm_mul_pd( t2, t3 );
        t2 = _mm_permute_pd( t2, 0b01 );
        t1 = _mm_fmaddsub_pd( t1, t3, t2 );
        mz = _mm_add_pd( mz, t1 );
      }

      t1 = _mm_loaddup_pd(&dd[i].re);
      t2 = _mm_loaddup_pd(&dd[i].im);
      t2 = _mm_mul_pd( t2, mz );
      t2 = _mm_permute_pd( t2, 0b01 );
      t1 = _mm_fmaddsub_pd( t1, mz, t2 );
      t1 = _mm_sub_pd( _mm_setzero_pd(), t1);  /* this line flips the sign of t1 */
      _mm_storeu_pd( &aa[6 * i + k].re, t1 );
    }
  }
}

void bck_house_avx512( complex_dble *aa, complex_dble *dd, double * rr )
{
  int i, j, k;
  complex_dble z;

  aa[35].re = dd[5].re;
  aa[35].im = dd[5].im;

  for (k = 4; k >= 0; k--) {
    z.re = dd[k].re;
    z.im = dd[k].im;
    dd[k].re = aa[6 * k + k].re;
    dd[k].im = aa[6 * k + k].im;
    aa[6 * k + k].re = z.re;
    aa[6 * k + k].im = z.im;

    for (j = (k + 1); j < 6; j++) {
      dd[j].re = aa[6 * j + k].re;
      dd[j].im = aa[6 * j + k].im;
      aa[6 * j + k].re = 0.0;
      aa[6 * j + k].im = 0.0;
    }

    for (i = 0; i < 6; i++) {
      __m128d mz, t1, t2, t3;
      mz = _mm_setzero_pd();

      for (j = k; j < 6; j++) {
        t1 = _mm_loaddup_pd(&aa[6 * i + j].re);
        t2 = _mm_loaddup_pd(&aa[6 * i + j].im);
        t3 = _mm_load_pd(&dd[j].re);
        t2 = _mm_mul_pd( t2, t3 );
        t2 = _mm_permute_pd( t2, 0b01 );
        t1 = _mm_fmaddsub_pd( t1, t3, t2 );
        mz = _mm_add_pd( mz, t1 );
      }

      t1 = _mm_loaddup_pd( rr+k );
      mz = _mm_mul_pd( mz, t1 );

      for (j = k; j < 6; j++) {
        t1 = _mm_loaddup_pd( &dd[j].re );
        t2 = _mm_loaddup_pd( &dd[j].im );
        t2 = _mm_mul_pd( t2, mz );
        t2 = _mm_permute_pd( t2, 0b01 );
        t1 = _mm_fmsubadd_pd( t1, mz, t2 );

        t2 = _mm_load_pd( &aa[6 * i + j].re );
        t1 = _mm_sub_pd( t2, t1 );
        _mm_storeu_pd( &aa[6 * i + j].re, t1 );
      }
    }
  }
}

#endif