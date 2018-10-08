/*******************************************************************************
*
* File salg_avx512.c
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementations of single precision linear algebra
* functions for spinors.
*
* See ../salg.c for more information and alternative
* implementations.
*
*******************************************************************************/

#ifdef AVX512

#include "global.h"
#include "linalg.h"
#include "mpi.h"
#include "sflds.h"

#include "avx512.h"

void mulc_spinor_add(int vol, spinor *s, spinor *r, complex z)
{
  spinor *sm;
  __m128 tr, ti;
  __m512 zr, zi, t1, t2;
  __m512 sign;

  sign = _mm512_set_ps( -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1 );
  sm = s + vol;

  tr = _mm_load_ps1( &z.re );
  ti = _mm_load_ps1( &z.im );
  zr = _mm512_broadcast_f32x4( tr );
  zi = _mm512_broadcast_f32x4( ti );

  zi = _mm512_mul_ps( zi, sign );

  for (; s < sm; s+=2) {
    t1 = _mm512_loadu_ps( &(*r).c1.c1.re );
    t2 = _mm512_mul_ps( zi, t1 );
    t2 = _mm512_permute_ps( t2, 0b10110001 );
    t2 = _mm512_fmadd_ps( zr, t1, t2 );
    t1 = _mm512_loadu_ps( &(*s).c1.c1.re );
    t1 = _mm512_add_ps( t1, t2 );
    _mm512_storeu_ps( &(*s).c1.c1.re, t1 );

    t1 = _mm512_loadu_ps( &(*r).c1.c1.re + 16 );
    t2 = _mm512_mul_ps( zi, t1 );
    t2 = _mm512_permute_ps( t2, 0b10110001 );
    t2 = _mm512_fmadd_ps( zr, t1, t2 );
    t1 = _mm512_loadu_ps( &(*s).c1.c1.re + 16 );
    t1 = _mm512_add_ps( t1, t2 );
    _mm512_storeu_ps( &(*s).c1.c1.re + 16, t1 );

    t1 = _mm512_loadu_ps( &(*r).c1.c1.re + 32 );
    t2 = _mm512_mul_ps( zi, t1 );
    t2 = _mm512_permute_ps( t2, 0b10110001 );
    t2 = _mm512_fmadd_ps( zr, t1, t2 );
    t1 = _mm512_loadu_ps( &(*s).c1.c1.re + 32 );
    t1 = _mm512_add_ps( t1, t2 );
    _mm512_storeu_ps( &(*s).c1.c1.re + 32, t1 );

    r += 2;
  }
}

#if __GNUC__ < 7
/* This function was implemented to gcc 7 */
extern __inline double _mm512_reduce_add_ps( __m512 a ) {
  float * d = (float *) &a;
  return d[0]+d[1]+d[2]+d[3]+d[4]+d[5]+d[6]+d[7]
        +d[8]+d[9]+d[10]+d[11]+d[12]+d[13]+d[14]+d[15] ;
}
#endif

complex spinor_prod(int vol, int icom, spinor *s, spinor *r )
{
  complex z;
  complex_dble v, w;
    spinor const *sm, *smb;
  __m512 tr, ti, s1, s2, s3, r1, r2, r3, sign;

  double x, y;

  x = 0.0;
  y = 0.0;
  sm = s + vol;


  while (s < sm) {
    smb = s + 8;
    if (smb > sm) {
      smb = sm;
    }

    tr = _mm512_setzero_ps();
    ti = _mm512_setzero_ps();

    for (; s < smb; s+=2) {
      s1 = _mm512_loadu_ps( &(*s).c1.c1.re );
      s2 = _mm512_loadu_ps( &(*s).c1.c1.re+16 );
      s3 = _mm512_loadu_ps( &(*s).c1.c1.re+32 );
      r1 = _mm512_loadu_ps( &(*r).c1.c1.re );
      r2 = _mm512_loadu_ps( &(*r).c1.c1.re+16 );
      r3 = _mm512_loadu_ps( &(*r).c1.c1.re+32 );

      tr = _mm512_fmadd_ps( s1, r1, tr );
      tr = _mm512_fmadd_ps( s2, r2, tr );
      tr = _mm512_fmadd_ps( s3, r3, tr );

      r1 = _mm512_permute_ps( r1, 0b10110001 );
      r2 = _mm512_permute_ps( r2, 0b10110001 );
      r3 = _mm512_permute_ps( r3, 0b10110001 );

      sign = _mm512_set_ps( -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1 );
      r1 = _mm512_mul_ps( r1, sign );
      r2 = _mm512_mul_ps( r2, sign );
      r3 = _mm512_mul_ps( r3, sign );

      ti = _mm512_fmadd_ps( s1, r1, ti );
      ti = _mm512_fmadd_ps( s2, r2, ti );
      ti = _mm512_fmadd_ps( s3, r3, ti );

      r += 2;
    }

    x += (double) _mm512_reduce_add_ps( tr );
    y += (double) _mm512_reduce_add_ps( ti );

  }

  v.re = x;
  v.im = y;

  if ((icom==1)&&(NPROC>1))
  {
     MPI_Reduce(&v.re,&w.re,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Bcast(&w.re,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
     z.re=(float)(w.re);
     z.im=(float)(w.im);
  }
  else
  {
     z.re=(float)(v.re);
     z.im=(float)(v.im);
  }
  return z;
}

#endif