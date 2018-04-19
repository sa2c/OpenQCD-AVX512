/*******************************************************************************
*
* File Dw_dble_avx512.c
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* AVX512 implementation of the O(a)-improved Wilson-Dirac operator D (double-
* precision programs).
*
* See ../Dw_dble.c for more details and alternative implementations
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
#include "dirac.h"
#include "global.h"

#define N0 (NPROC0 * L0)

typedef union
{
  spinor_dble s;
  weyl_dble w[2];
} spin_t;

#include "avx512.h"
#include "sse.h"

void doe_dble_avx512( const int *piup, const int *pidn, const su3_dble *u, const spinor_dble *pk, double coe, spin_t *rs)
{
  const spinor_dble *sp, *sm;
  const su3_dble *up;

  /* 512-bit wide stores for the spinor for each color */
  __m512d a1, a2, a3;
  __m512d b1, b2, b3;
  __m512d w1, w2, w3;
  __m512d t1, t2, t3, t4, t5, t6;

  __m128d tc;
  __m512d c512;

  /******************************* direction 0 *********************************/

  sp = pk + (*(piup++));
  sm = pk + (*(pidn++));

  _avx512_load_2_halfspinor_d( t1, t2, t3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_load_2_halfspinor_d( t4, t5, t6, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  a1 = _mm512_maskz_add_pd( 0b00001111, t1, t4 );
  a1 = _mm512_mask_sub_pd( a1, 0b11110000, t1, t4 );
  a2 = _mm512_maskz_add_pd( 0b00001111, t2, t5 );
  a2 = _mm512_mask_sub_pd( a2, 0b11110000, t2, t5 );
  a3 = _mm512_maskz_add_pd( 0b00001111, t3, t6 );
  a3 = _mm512_mask_sub_pd( a3, 0b11110000, t3, t6 );

  sp = pk + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pk + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );

  up = u;
  u += 1;
  avx512_su3_mul_quad_dble(  *up, *u, b1, b2, b3, a1, a2, a3  );

  _avx512_to_weyl_1( w1, b1 );
  _avx512_to_weyl_1( w2, b2 );
  _avx512_to_weyl_1( w3, b3 );


  /******************************* direction 1 *********************************/
  _avx512_load_2_halfspinor_d( t1, t2, t3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_load_2_halfspinor_d_reverse( t4, t5, t6, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  t4 = _mm512_permute_pd ( t4, 0b01010101 );
  a1 = _mm512_maskz_add_pd( 0b01011010, t1, t4 );
  a1 = _mm512_mask_sub_pd( a1, 0b10100101, t1, t4 );
  t5 = _mm512_permute_pd ( t5, 0b01010101 );
  a2 = _mm512_maskz_add_pd( 0b01011010, t2, t5 );
  a2 = _mm512_mask_sub_pd( a2, 0b10100101, t2, t5 );
  t6 = _mm512_permute_pd ( t6, 0b01010101 );
  a3 = _mm512_maskz_add_pd( 0b01011010, t3, t6 );
  a3 = _mm512_mask_sub_pd( a3, 0b10100101, t3, t6 );

  sp = pk + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pk + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );
  up = ++u;
  u += 1;

  avx512_su3_mul_quad_dble(  *up, *u, b1, b2, b3, a1, a2, a3  );

  _avx512_to_weyl_2( w1, b1 );
  _avx512_to_weyl_2( w2, b2 );
  _avx512_to_weyl_2( w3, b3 );

  /******************************* direction 2 *********************************/

  _avx512_load_2_halfspinor_d( t1, t2, t3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_load_2_halfspinor_d_reverse( t4, t5, t6, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  a1 = _mm512_maskz_add_pd( 0b11000011, t1, t4 );
  a1 = _mm512_mask_sub_pd( a1, 0b00111100, t1, t4 );
  a2 = _mm512_maskz_add_pd( 0b11000011, t2, t5 );
  a2 = _mm512_mask_sub_pd( a2, 0b00111100, t2, t5 );
  a3 = _mm512_maskz_add_pd( 0b11000011, t3, t6 );
  a3 = _mm512_mask_sub_pd( a3, 0b00111100, t3, t6 );

  sp = pk + (*(piup));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pk + (*(pidn));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );
  up = ++u;
  u += 1;

  avx512_su3_mul_quad_dble(  *up, *u, b1, b2, b3, a1, a2, a3  );

  _avx512_to_weyl_3( w1, b1 );
  _avx512_to_weyl_3( w2, b2 );
  _avx512_to_weyl_3( w3, b3 );


  /******************************* direction 3 *********************************/
  _avx512_load_2_halfspinor_d( t1, t2, t3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_load_2_halfspinor_d( t4, t5, t6, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  t4 = _mm512_permute_pd ( t4, 0b01010101 );
  a1 = _mm512_maskz_add_pd( 0b10010110, t1, t4 );
  a1 = _mm512_mask_sub_pd( a1, 0b01101001, t1, t4 );
  t5 = _mm512_permute_pd ( t5, 0b01010101 );
  a2 = _mm512_maskz_add_pd( 0b10010110, t2, t5 );
  a2 = _mm512_mask_sub_pd( a2, 0b01101001, t2, t5 );
  t6 = _mm512_permute_pd ( t6, 0b01010101 );
  a3 = _mm512_maskz_add_pd( 0b10010110, t3, t6 );
  a3 = _mm512_mask_sub_pd( a3, 0b01101001, t3, t6 );

  up = ++u;
  u += 1;
  avx512_su3_mul_quad_dble(  *up, *u, b1, b2, b3, a1, a2, a3  );

  _avx512_to_weyl_4( w1, b1 );
  _avx512_to_weyl_4( w2, b2 );
  _avx512_to_weyl_4( w3, b3 );

  tc =  _mm_load_sd( &coe );
  c512 = _mm512_broadcastsd_pd( tc );
  w1 = _mm512_mul_pd( c512, w1 );
  w2 = _mm512_mul_pd( c512, w2 );
  w3 = _mm512_mul_pd( c512, w3 );

  _avx512_store_2_halfspinor_d( w1, w2, w3, &rs->s.c1.c1.re, &rs->s.c3.c1.re );
}

void deo_dble_avx512( const int *piup, const int *pidn, const su3_dble *u,  spinor_dble *pl, double ceo, spin_t *rs)
{
  const su3_dble *up;
  spinor_dble *sp, *sm;

  /* 512-bit wide stores for the spinor for each color */
   __m512d a1, a2, a3;
   __m512d b1, b2, b3;
   __m512d w1, w2, w3;

   __m128d tc;
   __m512d c512;

  /******************************* direction 0 *********************************/

  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );

  _avx512_load_2_halfspinor_d( w1, w2, w3, &rs->s.c1.c1.re, &rs->s.c3.c1.re );

  tc =  _mm_load_sd( &ceo );
  c512 = _mm512_broadcastsd_pd( tc );
  w1 = _mm512_mul_pd( c512, w1 );
  w2 = _mm512_mul_pd( c512, w2 );
  w3 = _mm512_mul_pd( c512, w3 );

  _avx512_expand_weyl( a1, w1 )
  _avx512_expand_weyl( a2, w2 )
  _avx512_expand_weyl( a3, w3 )

  up = u;
  u += 1;
  avx512_su3_mul_quad_dble(  *u, *up, b1, b2, b3, a1, a2, a3  );

  _avx512_add_to_spinors( b1, b2, b3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_add_to_spinors_2( b1, b2, b3, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  /******************************* direction 1 *********************************/
  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );

  _avx512_expand_weyl_2( a1, w1 );
  _avx512_expand_weyl_2( a2, w2 );
  _avx512_expand_weyl_2( a3, w3 );

  up = ++u;
  u += 1;
  avx512_su3_mul_quad_dble(  *u, *up, b1, b2, b3, a1, a2, a3  );

  _avx512_add_to_spinors( b1, b2, b3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_add_to_spinors_3( b1, b2, b3, &(*sp).c3.c1.re, &(*sm).c3.c1.re );

  /******************************* direction 2 *********************************/
  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );

  _avx512_expand_weyl_3( a1, w1 );
  _avx512_expand_weyl_3( a2, w2 );
  _avx512_expand_weyl_3( a3, w3 );

  up = ++u;
  u += 1;
  avx512_su3_mul_quad_dble(  *u, *up, b1, b2, b3, a1, a2, a3  );

  _avx512_add_to_spinors( b1, b2, b3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_add_to_spinors_4( b1, b2, b3, &(*sp).c3.c1.re, &(*sm).c3.c1.re );


  /******************************* direction 3 *********************************/
  sp = pl + (*(piup++));
  _mm_prefetch( (char *) sp, _MM_HINT_T0 );
  sm = pl + (*(pidn++));
  _mm_prefetch( (char *) sm, _MM_HINT_T0 );

  _avx512_expand_weyl_4( a1, w1 );
  _avx512_expand_weyl_4( a2, w2 );
  _avx512_expand_weyl_4( a3, w3 );

  up = ++u;
  u += 1;
  avx512_su3_mul_quad_dble(  *u, *up, b1, b2, b3, a1, a2, a3  );

  _avx512_add_to_spinors( b1, b2, b3, &(*sp).c1.c1.re, &(*sm).c1.c1.re );
  _avx512_add_to_spinors_5( b1, b2, b3, &(*sp).c3.c1.re, &(*sm).c3.c1.re );
}
