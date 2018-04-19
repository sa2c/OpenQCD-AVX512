/*******************************************************************************
*
* File Dw_avx512.c
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementation of the O(a)-improved Wilson-Dirac operator D (single-
* precision programs).
*
* See ../Dw.c for more details and alternative implementations
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "utils.h"
#include "flags.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "sw_term.h"
#include "block.h"
#include "dirac.h"
#include "global.h"

#define N0 (NPROC0 * L0)

typedef union
{
  spinor s;
  weyl w[2];
} spin_t;

#include "avx512.h"
void doe_avx512(int *piup, int *pidn, su3 *u, spinor *pk, float coe, spin_t *rs)
{
  spinor *sp, *sm, *sp2, *sm2;
  su3 *up, *up1, *u1, *up2, *u2;

  /* 512-bit wide stores for the spinor for each color */
   __m512 a1, a2, a3;
   __m512 b1, b2, b3;
   __m256 c1, c2, c3;

   __m256 c256;

  /******************************* direction 0,1 **********************************/

  sp = pk + (*(piup++));
  sm = pk + (*(pidn++));
  sp2 = pk + (*(piup++));
  sm2 = pk + (*(pidn++));

  _avx512_load_4_halfspinor_f( a1, a2, a3,
                  &(*sp).c1.c1.re, &(*sm).c1.c1.re,
                  &(*sp2).c1.c1.re, &(*sm2).c1.c1.re );
  _avx512_load_4_halfspinor_f_reverse_up( b1, b2, b3,
                 &(*sp).c3.c1.re, &(*sm).c3.c1.re,
                 &(*sp2).c3.c1.re, &(*sm2).c3.c1.re );

  sp = pk + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pk + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );
  sp2 = pk + (*(piup));
  _mm_prefetch( (char *) sp2, _MM_HINT_T0 );
  sm2 = pk + (*(pidn));
  _mm_prefetch( (char *) sm2, _MM_HINT_T0 );

  up1 = u;
  u1  = u+1;
  up2 = u+2;
  u2  = u+3; u=u2;
  _avx512_dirac_combine_f_1( a1, b1 );
  _avx512_dirac_combine_f_1( a2, b2 );
  _avx512_dirac_combine_f_1( a3, b3 );

  avx512_su3_mixed_multiply_8( *up1, *u1, *up2, *u2, b1, b2, b3, a1, a2, a3 );

  _avx512_to_weyl_f_12( c1, b1 );
  _avx512_to_weyl_f_12( c2, b2 );
  _avx512_to_weyl_f_12( c3, b3 );

  /******************************* direction 2,3 *********************************/

  _avx512_load_4_halfspinor_f( a1, a2, a3,
                 &(*sp).c1.c1.re,&(*sm).c1.c1.re,
                 &(*sp2).c1.c1.re, &(*sm2).c1.c1.re );
  _avx512_load_4_halfspinor_f_reverse_dn( b1, b2, b3,
                 &(*sp).c3.c1.re, &(*sm).c3.c1.re,
                 &(*sp2).c3.c1.re, &(*sm2).c3.c1.re );

  _avx512_dirac_combine_f_2( a1, b1 );
  _avx512_dirac_combine_f_2( a2, b2 );
  _avx512_dirac_combine_f_2( a3, b3 );

  up1 = u+1;
  u1  = u+2;
  up2 = u+3;
  u2  = u+4;
  avx512_su3_mixed_multiply_8( *up1, *u1, *up2, *u2, b1, b2, b3, a1, a2, a3 );


  c256 = _mm256_broadcast_ss( &coe );

  _avx512_to_weyl_f_34( c1, b1 );
  _avx512_to_weyl_f_34( c2, b2 );
  _avx512_to_weyl_f_34( c3, b3 );

  c1 = _mm256_mul_ps( c1, c256);
  c2 = _mm256_mul_ps( c2, c256);
  c3 = _mm256_mul_ps( c3, c256);

  _avx512_write_6_hwv_f( c1, c2, c3,  &rs->s.c1.c1.re);
}

void deo_avx512(int *piup, int *pidn, su3 *u, spinor *pl, float ceo, spin_t *rs)
{
  spinor *sp, *sm, *sp2, *sm2;
  su3 *up, *up1, *u1, *up2, *u2;

  /* 512-bit wide stores for the spinor for each color */
  __m512 a1, a2, a3;
  __m512 b1, b2, b3;
  __m256 c1, c2, c3;

  __m256 c256;


  /******************************* direction 0 *********************************/

  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );
  sp2 = pl + (*(piup++));
  _mm_prefetch( (char *) sp2, _MM_HINT_T0 );
  sm2 = pl + (*(pidn++));
  _mm_prefetch( (char *) sm2, _MM_HINT_T0 );

  _avx512_load_6_hwv_f( c1,c2,c3, &rs->s.c1.c1.re );

  c256 = _mm256_broadcast_ss( &ceo );
  c1 = _mm256_mul_ps( c1, c256 );
  c2 = _mm256_mul_ps( c2, c256 );
  c3 = _mm256_mul_ps( c3, c256 );

  _avx512_to_dirac_f_1( a1, c1 );
  _avx512_to_dirac_f_1( a2, c2 );
  _avx512_to_dirac_f_1( a3, c3 );

  _avx512_to_dirac_f_2( a1, c1 );
  _avx512_to_dirac_f_2( a2, c2 );
  _avx512_to_dirac_f_2( a3, c3 );

  up1 = u;
  u1 = u+1;
  up2 = u+2;
  u2 = u+3;  u=u2;
  avx512_su3_mixed_multiply_8( *u1, *up1, *u2, *up2, b1, b2, b3, a1, a2, a3 );

  _avx512_load_4_halfspinor_f( a1, a2, a3,
                 &(*sm).c1.c1.re, &(*sp).c1.c1.re,
                 &(*sm2).c1.c1.re, &(*sp2).c1.c1.re );
  a1 = _mm512_add_ps( a1, b1 );
  a2 = _mm512_add_ps( a2, b2 );
  a3 = _mm512_add_ps( a3, b3 );
  _avx512_write_4_halfspinor_f( a1, a2, a3,
                  &(*sm).c1.c1.re, &(*sp).c1.c1.re,
                  &(*sm2).c1.c1.re, &(*sp2).c1.c1.re );

  _avx512_load_4_halfspinor_f_reverse_up( a1, a2, a3,
                 &(*sm).c3.c1.re, &(*sp).c3.c1.re,
                 &(*sm2).c3.c1.re, &(*sp2).c3.c1.re );
  _avx512_dirac_combine_f_3( a1, b1 );
  _avx512_dirac_combine_f_3( a2, b2 );
  _avx512_dirac_combine_f_3( a3, b3 );
  _avx512_write_4_halfspinor_f_reverse_up( a1, a2, a3,
                  &(*sm).c3.c1.re, &(*sp).c3.c1.re,
                  &(*sm2).c3.c1.re, &(*sp2).c3.c1.re );

  /******************************* direction 2 *********************************/

  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );
  sp2 = pl + (*(piup++));
  _mm_prefetch( (char *) sp2, _MM_HINT_T0 );
  sm2 = pl + (*(pidn++));
  _mm_prefetch( (char *) sm2, _MM_HINT_T0 );

  _avx512_to_dirac_f_3( a1, c1 );
  _avx512_to_dirac_f_3( a2, c2 );
  _avx512_to_dirac_f_3( a3, c3 );

  up1 = u+1;
  u1 = u+2;
  up2 = u+3;
  u2 = u+4;
  avx512_su3_mixed_multiply_8( *u1, *up1, *u2, *up2, b1, b2, b3, a1, a2, a3 );

  _avx512_load_4_halfspinor_f( a1, a2, a3, &(*sm).c1.c1.re, &(*sp).c1.c1.re, &(*sm2).c1.c1.re, &(*sp2).c1.c1.re );
  a1 = _mm512_add_ps( a1, b1 );
  a2 = _mm512_add_ps( a2, b2 );
  a3 = _mm512_add_ps( a3, b3 );
  _avx512_write_4_halfspinor_f( a1, a2, a3, &(*sm).c1.c1.re, &(*sp).c1.c1.re, &(*sm2).c1.c1.re, &(*sp2).c1.c1.re );

  _avx512_load_4_halfspinor_f_reverse_dn( a1, a2, a3, &(*sm).c3.c1.re, &(*sp).c3.c1.re, &(*sm2).c3.c1.re, &(*sp2).c3.c1.re );
  _avx512_dirac_combine_f_4( a1, b1 );
  _avx512_dirac_combine_f_4( a2, b2 );
  _avx512_dirac_combine_f_4( a3, b3 );
  _avx512_write_4_halfspinor_f_reverse_dn( a1, a2, a3, &(*sm).c3.c1.re, &(*sp).c3.c1.re, &(*sm2).c3.c1.re, &(*sp2).c3.c1.re );
}
