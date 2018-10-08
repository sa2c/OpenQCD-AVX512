/*******************************************************************************
*
* File salg_dble_avx512.c
* 
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementations of single precision linear algebra
* functions for spinors.
*
* See ../salg_dble.c for more information and alternative
* implementations.
*
*******************************************************************************/

#ifdef AVX512

#include "global.h"
#include "linalg.h"
#include "mpi.h"
#include "sflds.h"

#include "avx512.h"
#if __GNUC__ < 7
/* This function was implemented to gcc 7 */
extern __inline double _mm512_reduce_add_pd( __m512d a ) {
  double * d = (double *) &a;
  return d[0]+d[1]+d[2]+d[3]+d[4]+d[5]+d[6]+d[7] ;
}
#endif

complex_dble spinor_prod_dble_avx512( spinor_dble const *s, spinor_dble const *smb, spinor_dble const *r)
{
  __m512d tr, ti, s1, s2, s3, r1, r2, r3, sign;
  sign = _mm512_set_pd( -1,1,-1,1,-1,1,-1,1 );
  complex_dble z;

  tr = _mm512_setzero_pd();
  ti = _mm512_setzero_pd();
  for (; s < smb; s++) {
    s1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    s2 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    s3 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    r1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    r2 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    r3 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    tr = _mm512_fmadd_pd( s1, r1, tr );
    tr = _mm512_fmadd_pd( s2, r2, tr );
    tr = _mm512_fmadd_pd( s3, r3, tr );
    r1 = _mm512_permute_pd( r1, 0b01010101 );
    r2 = _mm512_permute_pd( r2, 0b01010101 );
    r3 = _mm512_permute_pd( r3, 0b01010101 );
    sign = _mm512_set_pd( -1,1,-1,1,-1,1,-1,1 );
    r1 = _mm512_mul_pd( r1, sign );
    r2 = _mm512_mul_pd( r2, sign );
    r3 = _mm512_mul_pd( r3, sign );
    ti = _mm512_fmadd_pd( s1, r1, ti );
    ti = _mm512_fmadd_pd( s2, r2, ti );
    ti = _mm512_fmadd_pd( s3, r3, ti );
    r += 1;
  }
  z.re = _mm512_reduce_add_pd( tr );
  z.im = _mm512_reduce_add_pd( ti );
  return z;
}

double spinor_prod_re_dble_avx512( spinor_dble const *s, spinor_dble const *smb, spinor_dble const *r)
{
  __m512d tr, ti, s1, s2, s3, r1, r2, r3;
  float c;
  tr = _mm512_setzero_pd();
  ti = _mm512_setzero_pd();
  for (; s < smb; s++) {
    s1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    s2 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    s3 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    r1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    r2 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    r3 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    tr = _mm512_fmadd_pd( s1, r1, tr );
    tr = _mm512_fmadd_pd( s2, r2, tr );
    tr = _mm512_fmadd_pd( s3, r3, tr );
    r1 = _mm512_permute_pd( r1, 0b01010101 );
    r2 = _mm512_permute_pd( r2, 0b01010101 );
    r3 = _mm512_permute_pd( r3, 0b01010101 );
    r += 1;
  }
  c = _mm512_reduce_add_pd( tr );
  return c;
}

complex_dble spinor_prod5_dble_avx512(spinor_dble const *s, spinor_dble const *smb, spinor_dble const *r)
{
  __m512d tr, ti, s1, s2, s3, r1, r2, r3, sign;
  complex_dble z;

  tr = _mm512_setzero_pd();
  ti = _mm512_setzero_pd();
  for (; s < smb; s++) {
    s1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    s2 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    s3 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    r1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    r2 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    r3 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    sign = _mm512_set_pd( -1,-1,-1,-1,1,1,1,1 );
    s2 = _mm512_mul_pd( s2, sign );
    tr = _mm512_fmadd_pd( s1, r1, tr );
    tr = _mm512_fmadd_pd( s2, r2, tr );
    tr = _mm512_fnmadd_pd( s3, r3, tr );
    r1 = _mm512_permute_pd( r1, 0b01010101 );
    r2 = _mm512_permute_pd( r2, 0b01010101 );
    r3 = _mm512_permute_pd( r3, 0b01010101 );
    sign = _mm512_set_pd( -1,1,-1,1,-1,1,-1,1 );
    r1 = _mm512_mul_pd( r1, sign );
    r2 = _mm512_mul_pd( r2, sign );
    r3 = _mm512_mul_pd( r3, sign );
    ti = _mm512_fmadd_pd( s1, r1, ti );
    ti = _mm512_fmadd_pd( s2, r2, ti );
    ti = _mm512_fnmadd_pd( s3, r3, ti );
    r += 1;
  }
  z.re = _mm512_reduce_add_pd( tr );
  z.im = _mm512_reduce_add_pd( ti );
  return z;
}

double norm_square_dble_avx512(spinor_dble const *s, spinor_dble const *smb)
{
  __m512d tmp, s1, s2, s3;
  tmp = _mm512_setzero_pd();
  for (; s < smb; s++) {
    s1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    s2 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    s3 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    tmp = _mm512_fmadd_pd( s1, s1, tmp );
    tmp = _mm512_fmadd_pd( s2, s2, tmp );
    tmp = _mm512_fmadd_pd( s3, s3, tmp );
  }
  return _mm512_reduce_add_pd( tmp );
}



void mulc_spinor_add_dble(int vol, spinor_dble *s, spinor_dble *r,
                          complex_dble z)
{
  spinor_dble *sm;
  __m128d tr, ti;
  __m512d zr, zi, t1, t2;

  tr = _mm_load_sd( &z.re );
  ti = _mm_load_sd( &z.im );
  zr = _mm512_broadcastsd_pd( tr );
  zi = _mm512_broadcastsd_pd( ti );

  sm = s + vol;

  for (; s < sm; s++) {
    t1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    t2 = _mm512_fmaddsub_pd( zr, t1, t2 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    t2 = _mm512_fmaddsub_pd( zr, t1, t2 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+8, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    t2 = _mm512_fmaddsub_pd( zr, t1, t2 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+16, t1 );

    r += 1;
  }
}


void mulr_spinor_add_dble(int vol, spinor_dble *s, spinor_dble *r,
                          double c)
{
  spinor_dble *sm;
  __m128d t128;
  __m512d tc, t1, t2;

  t128 = _mm_load_sd( &c );
  tc = _mm512_broadcastsd_pd( t128 );

  sm = s + vol;

  for (; s < sm; s++) {
    t1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    t2 = _mm512_mul_pd( tc, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    t2 = _mm512_mul_pd( tc, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+8, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    t2 = _mm512_mul_pd( tc, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    t1 = _mm512_add_pd( t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+16, t1 );

    r += 1;
  }
}


void combine_spinor_dble(int vol, spinor_dble *s, spinor_dble *r,
                         double cs, double cr)
{
  spinor_dble *sm;
  __m128d ts128, tr128;
  __m512d tcs, tcr, t1, t2;

  ts128 = _mm_load_sd( &cs );
  tr128 = _mm_load_sd( &cr );
  tcs = _mm512_broadcastsd_pd( ts128 );
  tcr = _mm512_broadcastsd_pd( tr128 );

  sm = s + vol;

  for (; s < sm; s++) {
    t1 = _mm512_loadu_pd( &(*r).c1.c1.re );
    t2 = _mm512_mul_pd( tcr, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    t1 = _mm512_fmadd_pd( tcs, t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+8 );
    t2 = _mm512_mul_pd( tcr, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    t1 = _mm512_fmadd_pd( tcs, t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+8, t1 );

    t1 = _mm512_loadu_pd( &(*r).c1.c1.re+16 );
    t2 = _mm512_mul_pd( tcr, t1 );
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    t1 = _mm512_fmadd_pd( tcs, t1, t2 );
    _mm512_storeu_pd( &(*s).c1.c1.re+16, t1 );

    r += 1;
  }
}


void scale_dble(int vol, double c, spinor_dble *s)
{
  spinor_dble *sm;
  __m128d t128;
  __m512d tc, t1;

  t128 = _mm_load_sd( &c );
  tc = _mm512_broadcastsd_pd( t128 );

  sm = s + vol;

  for (; s < sm; s++) {
    t1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    t1 = _mm512_mul_pd( tc, t1 );
    _mm512_storeu_pd( &(*s).c1.c1.re, t1 );

    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+8 );
    t1 = _mm512_mul_pd( tc, t1 );
    _mm512_storeu_pd( &(*s).c1.c1.re+8, t1 );

    t1 = _mm512_loadu_pd( &(*s).c1.c1.re+16 );
    t1 = _mm512_mul_pd( tc, t1 );
    _mm512_storeu_pd( &(*s).c1.c1.re+16, t1 );
  }
}

void rotate_dble_avx512(int n, int ix, spinor_dble **ppk, spinor_dble *psi, complex_dble const *v)
{
  spinor_dble *pk, *pj;
  complex_dble const *z;
  int k,j;

  for (k = 0; k < n; k++) {
    __m128d tr, ti;
    __m512d zr,zi, t1,t2, p1,p2,p3, sign;

    pk = psi + k;
    pj = ppk[0] + ix;
    z = v + k;

    tr = _mm_load_pd1( &z->re );
    ti = _mm_load_pd1( &z->im );
    zr = _mm512_broadcastsd_pd( tr );
    zi = _mm512_broadcastsd_pd( ti );

    sign = _mm512_set_pd( -1,1,-1,1,-1,1,-1,1 );
    zi = _mm512_mul_pd( zi, sign );

    t1 = _mm512_loadu_pd( &(*pj).c1.c1.re );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    p1 = _mm512_fmadd_pd( zr, t1, t2 );

    t1 = _mm512_loadu_pd( &(*pj).c1.c1.re+8 );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    p2 = _mm512_fmadd_pd( zr, t1, t2 );

    t1 = _mm512_loadu_pd( &(*pj).c1.c1.re+16 );
    t2 = _mm512_mul_pd( zi, t1 );
    t2 = _mm512_permute_pd( t2, 0b01010101 );
    p3 = _mm512_fmadd_pd( zr, t1, t2 );

    for (j = 1; j < n; j++) {
      pj = ppk[j] + ix;
      z += n;

      tr = _mm_load_pd1( &z->re );
      ti = _mm_load_pd1( &z->im );
      zr = _mm512_broadcastsd_pd( tr );
      zi = _mm512_broadcastsd_pd( ti );
      zi = _mm512_mul_pd( zi, sign );

      t1 = _mm512_loadu_pd( &(*pj).c1.c1.re );
      t2 = _mm512_mul_pd( zi, t1 );
      t2 = _mm512_permute_pd( t2, 0b01010101 );
      t1 = _mm512_fmadd_pd( zr, t1, t2 );
      p1 = _mm512_add_pd( p1, t1 );

      t1 = _mm512_loadu_pd( &(*pj).c1.c1.re+8 );
      t2 = _mm512_mul_pd( zi, t1 );
      t2 = _mm512_permute_pd( t2, 0b01010101 );
      t1 = _mm512_fmadd_pd( zr, t1, t2 );
      p2 = _mm512_add_pd( p2, t1 );

      t1 = _mm512_loadu_pd( &(*pj).c1.c1.re+16 );
      t2 = _mm512_mul_pd( zi, t1 );
      t2 = _mm512_permute_pd( t2, 0b01010101 );
      t1 = _mm512_fmadd_pd( zr, t1, t2 );
      p3 = _mm512_add_pd( p3, t1 );
    }

    _mm512_storeu_pd( &(*pk).c1.c1.re, p1 );
    _mm512_storeu_pd( &(*pk).c1.c1.re+8, p2 );
    _mm512_storeu_pd( &(*pk).c1.c1.re+16, p3 );
  }
}


void mulg5_dble(int vol, spinor_dble *s)
{
  spinor_dble *sm;

  sm = s + vol;

  for (; s < sm; s++) {
    __m512d s1;
    __m256d s2;

    s1 = _mm512_loadu_pd( &(*s).c1.c1.re+12 );
    s1 = _mm512_sub_pd( _mm512_setzero_pd(), s1 );
    _mm512_storeu_pd( &(*s).c1.c1.re+12, s1 );

    s2 = _mm256_loadu_pd( &(*s).c1.c1.re+20 );
    s2 = _mm256_sub_pd( _mm256_setzero_pd(), s2 );
    _mm256_storeu_pd( &(*s).c1.c1.re+20, s2 );
  }
}

void mulmg5_dble(int vol, spinor_dble *s)
{
  spinor_dble *sm;

  sm = s + vol;

  for (; s < sm; s++) {
    __m512d s1;
    __m256d s2;

    s1 = _mm512_loadu_pd( &(*s).c1.c1.re );
    s1 = _mm512_sub_pd( _mm512_setzero_pd(), s1 );
    _mm512_storeu_pd( &(*s).c1.c1.re, s1 );

    s2 = _mm256_loadu_pd( &(*s).c1.c1.re+8 );
    s2 = _mm256_sub_pd( _mm256_setzero_pd(), s2 );
    _mm256_storeu_pd( &(*s).c1.c1.re+8, s2 );
  }
}

#endif